import os
import subprocess
import torch
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import List, Dict, Tuple, Optional

def one_hot_encode_sequence(sequence: str, max_length: int = 50) -> np.ndarray:
    """将DNA序列转换为独热编码"""
    nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 3}
    encoded = np.zeros((4, max_length), dtype=np.float32)
    
    # 截断或填充到指定长度
    sequence = sequence[:max_length].upper()
    for i, nucleotide in enumerate(sequence):
        if nucleotide in nucleotide_to_index:
            encoded[nucleotide_to_index[nucleotide], i] = 1.0
    
    return encoded

def run_command(cmd: List[str], cwd: Optional[str] = None) -> Tuple[str, str, int]:
    """运行系统命令并返回输出"""
    try:
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            cwd=cwd,
            check=False
        )
        return result.stdout, result.stderr, result.returncode
    except Exception as e:
        return "", str(e), 1

def parse_minced_output(minced_gff_file: str) -> List[Dict]:
    """解析MinCED的GFF输出文件"""
    crispr_arrays = []
    
    if not os.path.exists(minced_gff_file):
        return crispr_arrays
    
    with open(minced_gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) >= 9 and fields[2] == 'repeat_region':
                array_info = {
                    'seq_id': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'score': float(fields[5]) if fields[5] != '.' else 0.0,
                    'strand': fields[6],
                    'attributes': fields[8]
                }
                crispr_arrays.append(array_info)
    
    return crispr_arrays

def extract_repeat_sequences(fasta_file: str, crispr_arrays: List[Dict]) -> List[Dict]:
    """从FASTA文件中提取重复序列"""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    
    repeat_data = []
    for array in crispr_arrays:
        seq_id = array['seq_id']
        if seq_id in sequences:
            # 提取重复序列（简化：取第一个重复单元）
            start = array['start'] - 1  # 转换为0-based索引
            end = min(start + 50, len(sequences[seq_id]))  # 最多取50bp
            repeat_seq = sequences[seq_id][start:end]
            
            repeat_data.append({
                'array_id': f"{seq_id}_{array['start']}_{array['end']}",
                'sequence': repeat_seq,
                'seq_id': seq_id,
                'start': array['start'],
                'end': array['end']
            })
    
    return repeat_data

def load_fasta_sequences(fasta_file: str) -> Dict[str, str]:
    """加载FASTA文件中的所有序列"""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def save_results_to_csv(results: List[Dict], output_file: str):
    """保存结果到CSV文件"""
    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")