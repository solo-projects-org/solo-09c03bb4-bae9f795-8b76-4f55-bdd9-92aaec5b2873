# FASTQ 混样拆分与建库接头清洗服务

基于 FastAPI 的网页服务，分两阶段处理测序数据：

1. **混样拆分（demultiplex）**：把混样的 R1 / R2 / I1 / I2 四行制
   `FASTQ.gz` 按样本索引拆回各样本，产出可下载的 ZIP；
2. **建库接头清洗（cutadapt trimming）**：对已完成批次的**非空样本**，
   用 cutadapt 去掉 R1、R2 的 3′ 建库接头，任一端过短则整对丢弃，产出
   全新的双端文件与清洗 ZIP。**原拆分文件始终只读保留。**

## 一、拆分规则（接口与规则保持不变）

- 上传文件：`R1`、`R2`、`I1`、`I2` 四个 gzip 压缩的四行制 FASTQ，
  以及三列样本表 CSV（表头必须为 `sample_id,i7,i5`）。
- 样本索引等长、`i7+i5` 组合唯一；`i5` 按**匹配方向**填写（服务端不再做反向互补）。
- 四路读段名取**首个空白符之前**的部分，并去掉末尾 `/1`、`/2` 后对齐；
  任一条目四路名字不一致即整批失败。
- `I1` 比对样本 `i7`、`I2` 比对样本 `i5`，每端汉明距离 **≤ 1** 才进入候选；
  在候选中取**总距离最小且唯一**的样本。无候选或最优并列时，该读段进入
  `undetermined`（未分配）组。
- 以下情况**整批失败**（状态 `failed`，不提供下载）：四路条数不同、
  读段名错序 / 对不齐、序列与质量不等长、FASTQ 结构损坏、gzip 流损坏、样本表不合法。
- 逐条流式处理；只有**全部读段校验并保存成功**后才生成并开放拆分 ZIP。

## 二、接头清洗规则（cutadapt）

- **可选对象**：仅已完成（`completed`）拆分批次中**非空**的样本
  （含 `undetermined`，只要它有读段）。网页下拉框与 `GET /api/trimming/targets`
  只列出这些样本。
- **输入参数**：R1 的 3′ 接头（cutadapt `-a`）、R2 的 3′ 接头（`-A`）、
  最短长度（`-m`）。接头允许 IUPAC 核苷酸字符（A/C/G/T/N/…），
  最短长度为 ≥1 的整数；非法参数返回 `400`。
- **实际命令**（默认匹配规则，不改变错误率 / 重叠默认值）：

  ```bash
  cutadapt -j 2 \
    -a <adapter_r1> -A <adapter_r2> \
    --pair-filter=any -m <min_length> \
    -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
    --report=full <source_R1> <source_R2>
  ```

- `-j 2` 双线程；`--pair-filter=any`：**任一端**裁剪后短于最短长度，
  **整对丢弃**（两端同步过滤，保证双端仍然一一配对）。
- cutadapt 完成后服务端再校验：两个输出均存在且为合法 gzip FASTQ、
  R1/R2 记录条数一致；并校验生成的 ZIP 完整性。全部成功才发布。
- **原拆分文件只读，永不修改**；清洗产物写到独立的暂存目录与 ZIP。

### 重提、取代与取消（同样本）

- **同样本重提参数一律受理**（返回 `201`），并**取代**该样本运行中或排队中的清洗。
  每个样本只有一条串行 worker 线程：旧尝试的 cutadapt 进程组被终止、
  其整个临时目录被删除，**确认旧尝试的计算进程全部退出、临时文件清完后**，
  worker 才启动最后受理的参数。
- 被取代（`superseded`）或被取消（`cancelled`）的尝试**绝不发布**。
- **包的保留与原子替换**：新双端文件与 ZIP **完整成功前**，旧清洗包一直保留；
  新尝试失败或取消也保留旧包。发布采用指针原子切换，切换后才删除旧包目录。
  页面会对失败 / 取消 / 被取代且旧包仍在的尝试，注明**旧包参数**与
  **“新参数未生效”**。
- **下载一致性**：下载响应持有该版本 ZIP 的文件句柄直到传完；即使期间发生
  换包（旧目录被删），已开始的下载仍能完整收到**同一版**文件。
- **取消语义**：取消接口会阻塞到 cutadapt 进程组真正退出、临时目录清理完成后
  才确认。该路由为同步函数，由 FastAPI 线程池执行，等待**不会卡住事件循环、
  其他样本或原拆分服务**。不同样本的清洗在各自独立线程运行，互不阻塞。
- 清洗状态**只存内存**，进程重启即清空，启动时同时清空清洗暂存目录。

## 目录结构

```
app/
  main.py        FastAPI 路由（拆分 + 清洗全部 HTTP 接口、首页）
  jobs.py        拆分内存态作业注册表、双批信号量、取消与清理
  trim_jobs.py   清洗注册表：每样本串行 worker、取代/取消、原子换包、下载句柄
  trimming.py    cutadapt 调用（-j 2、3′ 接头、pair-filter=any）、校验、ZIP、报告解析
  pipeline.py    流式拆分主流程、分组 gzip 输出、ZIP 打包
  fastq_io.py    gzip FASTQ 流式读取器与损坏检测
  demux.py       样本表解析、汉明距离、候选/并列判定
  errors.py      领域异常（含 TrimmingError）
  templates/     单页前端（拆分 + 清洗）
scripts/
  generate_samples.py   生成四套演示数据
samples/
  normal/      正常拆分样例（全部唯一命中，含 1 错配）
  ambiguous/   歧义样例（含无候选与并列最优，各 1 条进 undetermined）
  corrupt/     末尾截断样例（R2 gzip 尾部损坏，整批失败）
  trimming/    清洗样例：R1/R2 含已知 3′ 接头、内部接头、过短读段与未分配读段
tests/         pytest 测试（50 项）
```

## 本地运行（无 Docker）

```bash
python3 -m pip install -r requirements.txt   # 含 cutadapt
uvicorn app.main:app --host 0.0.0.0 --port 8000
# 可选环境变量：
#   export DEMUX_JOBS_ROOT=/data/jobs        # 拆分暂存目录，默认 /data/jobs
#   export DEMUX_TRIM_ROOT=/data/trimming    # 清洗暂存目录，默认 /data/trimming
```

浏览器打开 <http://localhost:8000>：

1. 从 `samples/<场景>/` 选择五个文件上传，完成拆分；
2. 在“建库接头清洗”卡片的下拉框选择已完成批次的非空样本；
3. 清洗样例 `samples/trimming/` 的推荐参数（见该目录 `TRIMMING_README.txt`）：
   - R1 3′ 接头：`AGATCGGAAGAG`
   - R2 3′ 接头：`CTGTCTCTTATA`
   - 最短长度：`10`
   - 预期：SEQ-0003 一对因裁剪后仅 8 bp 被整对过滤，其余读段保留并去接头。

重新生成演示数据：

```bash
python3 -I scripts/generate_samples.py samples
```

## Docker Compose 运行

```bash
docker compose up --build
# 浏览器打开 http://localhost:8000
```

容器以非 root 用户运行；cutadapt 通过 pip wheel 安装（含其原生加速扩展）。
拆分与清洗分别写入命名卷 `demux-jobs`（`/data/jobs`）与
`demux-trimming`（`/data/trimming`），每次启动自动清空；容器内置健康检查。

## HTTP 接口

### 拆分（原有，行为不变）

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `GET`  | `/` | 上传、拆分监控与清洗网页 |
| `GET`  | `/api/health` | 健康检查 |
| `POST` | `/api/jobs` | multipart 上传 5 个文件并启动拆分，`201`；并发已满 `409` |
| `GET`  | `/api/jobs` | 全部拆分作业快照 |
| `GET`  | `/api/jobs/{id}` | 单个拆分作业快照 |
| `POST` | `/api/jobs/{id}/cancel` | 取消拆分并等待实际停止后返回 |
| `GET`  | `/api/jobs/{id}/download` | 仅完成态可下载拆分 ZIP，否则 `404` |

### 接头清洗（新增）

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `GET`  | `/api/trimming/targets` | 列出已完成批次的非空样本及其清洗状态 |
| `POST` | `/api/trimming` | 表单参数 `job_id,sample_id,adapter_r1,adapter_r2,min_length`；同样本重提即取代，`201`；参数非法 `400`，目标不存在/为空 `404` |
| `GET`  | `/api/trimming` | 全部清洗尝试快照（含当前发布包参数 `active_published`） |
| `GET`  | `/api/trimming/{id}` | 单个清洗尝试快照（状态、参数、cutadapt 统计、错误、是否发布） |
| `POST` | `/api/trimming/{id}/cancel` | 取消并等待 cutadapt 进程退出与清理完成后返回 |
| `GET`  | `/api/trimming/{id}/download` | 仅**已发布**的成功尝试可下载清洗 ZIP，否则 `404`；换包不影响进行中的下载 |

清洗尝试状态：`queued`（排队）、`running`（cutadapt 运行中）、
`completed`（成功并发布）、`failed`（失败，旧包保留）、
`cancelled`（用户取消，旧包保留）、`superseded`（被更新参数取代，不发布）。

清洗 ZIP 内容：

```
trimmed_R1.fastq.gz     # 去接头后的双端 R1
trimmed_R2.fastq.gz     # 去接头后的双端 R2（与 R1 配对、条数一致）
trim_params.txt         # 本次参数、-j 2、pair-filter=any 与完整命令行
cutadapt_report.txt     # cutadapt full 报告（含处理/保留/过滤/命中接头统计）
```

## 测试

```bash
python3 -m pytest tests/ -q
```

覆盖：原有索引匹配 / 样本表 / gzip·FASTQ 损坏 / 双批上限 / 取消名额释放等
27 项；清洗新增 23 项——参数校验、真实 cutadapt 去接头与短读段整对过滤、
undetermined 可清洗、原拆分文件不变、重提取代且旧尝试不发布、
取消等待进程退出并保留旧包、失败保留旧包且回显旧参数、取消不阻塞其他样本、
换包不截断进行中的下载、临时目录清理、`-j 2`/`-a`/`-A`/`--pair-filter=any`
命令行断言与进程组终止信号。
