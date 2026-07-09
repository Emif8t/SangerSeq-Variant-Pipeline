from config import AB1_FOLDER

from scripts.preprocessing import *

abi_records = load_ab1_files(AB1_FOLDER)

processed_reads = prepare_reads(abi_records)

qc_summary = calculate_qc_metrics(processed_reads)

save_qc_summary(qc_summary)
