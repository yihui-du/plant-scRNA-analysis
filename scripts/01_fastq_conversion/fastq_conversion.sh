#!/usr/bin/env bash
set -euo pipefail

# ====== 你需要改的参数 ======
ACC_LIST="accessions.txt"     # SRR/ERR/DRR 列表
OUTDIR="./fastq_raw"
THREADS=8                      # 下载/转换线程
TMPDIR="./_sra_tmp"
# ===========================

mkdir -p "${OUTDIR}" "${TMPDIR}"

# 需要：prefetch、fasterq-dump（SRA Toolkit）
# 建议先：vdb-config -i  (把缓存/下载目录配置好)

while read -r ACC; do
  [[ -z "${ACC}" ]] && continue
  echo "==> Processing ${ACC}"

  # 1) 下载 .sra 到本地缓存（或指定目录）
  prefetch "${ACC}" --max-size 200G --output-directory "${TMPDIR}"

  # 2) 转 FASTQ（--split-files 生成 _1/_2；单端会只出一个）
  fasterq-dump "${TMPDIR}/${ACC}/${ACC}.sra" \
    --split-files \
    --threads "${THREADS}" \
    --outdir "${OUTDIR}" \
    --temp "${TMPDIR}"

  # 3) 压缩（10x/cellranger 通常吃 .fastq.gz）
  gzip -f "${OUTDIR}/${ACC}"*.fastq

done < "${ACC_LIST}"

echo "Done. FASTQ in ${OUTDIR}"
