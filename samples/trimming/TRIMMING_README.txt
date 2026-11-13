Adapter trimming demo parameters (cutadapt -j 2, default matching):
  adapter_r1 = AGATCGGAAGAG
  adapter_r2 = CTGTCTCTTATA
  min_length = 10
Expected: SEQ-0003 pair is filtered (8 bp insert on both ends);
all other pairs survive with their 3' adapters removed.
