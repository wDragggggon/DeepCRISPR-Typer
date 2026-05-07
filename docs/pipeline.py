import os
import tempfile
from typing import List, Dict, Optional
from pathlib import Path
import concurrent.futures # 在文件顶部加上这个

from config.paths_config import MINCED_PATH, PRODIGAL_PATH, INPUT_DIR, OUTPUT_DIR
from src.utils import run_command, parse_minced_output, extract_repeat_sequences, load_fasta_sequences, save_results_to_csv
from src.repeat_classifier import RepeatClassifier
from src.cas_protein_classifier import CasProteinClassifier
from src.hmm_scanner import HMMScanner
from src.rule_based_classifier import RuleBasedClassifier

# 【新增】：必须放在类外面的顶层，否则多进程无法 Pickle！
def _run_prodigal_chunk_worker(chunk_fna: Path) -> Path:
    import subprocess
    import os
    from config.paths_config import PRODIGAL_PATH  # 确保能拿到路径
    
    # 优雅地把 .fna 后缀换成 .faa
    chunk_faa = chunk_fna.with_suffix('.faa') 
    cmd = [
        PRODIGAL_PATH,
        "-i", str(chunk_fna),
        "-a", str(chunk_faa),
        "-p", "meta",
        "-q"
    ]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return chunk_faa

class DeepCRISPRTyperPipeline:
    """DeepCRISPR-Typer主pipeline"""
    
    def __init__(self):
        self.repeat_classifier = RepeatClassifier()
        self.cas_classifier = CasProteinClassifier()
        self.hmm_scanner = HMMScanner()
        self.rule_classifier = RuleBasedClassifier()
    
    def run(self, input_fasta: str, output_prefix: str = None) -> List[Dict]:
        """运行完整的CRISPR-Cas分类pipeline"""
        
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
        
        if output_prefix is None:
            output_prefix = Path(input_fasta).stem
        
        output_dir = OUTPUT_DIR / output_prefix
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"🚀 Starting DeepCRISPR-Typer pipeline for {input_fasta}")
        print(f"📁 Output directory: {output_dir}")
        
        
        # # Step 1: CRISPR阵列检测 (MinCED)
        # print("\n🔍 Step 1: Detecting CRISPR arrays with MinCED...")
        # crispr_arrays = self._detect_crispr_arrays(input_fasta, output_dir)
        # if not crispr_arrays:
        #     print("⚠️  No CRISPR arrays detected!")
        #     return []
        
        # # Step 2: 提取重复序列
        # print("\n🧬 Step 2: Extracting repeat sequences...")
        # repeat_data = extract_repeat_sequences(input_fasta, crispr_arrays)
        # if not repeat_data:
        #     print("⚠️  No repeat sequences extracted!")
        #     return []
        
        # # Step 3: 重复序列分类
        # print("\n🧠 Step 3: Classifying repeat sequences...")
        # repeat_sequences = [r['sequence'] for r in repeat_data]
        # repeat_predictions = self.repeat_classifier.predict(repeat_sequences)
        
        # # Step 4: Cas蛋白预测 (Prodigal)
        # print("\n🔬 Step 4: Predicting protein coding genes with Prodigal...")
        # protein_fasta = self._predict_proteins(input_fasta, output_dir)
        # if not protein_fasta or not os.path.exists(protein_fasta):
        #     print("⚠️  Protein prediction failed!")
        #     protein_fasta = None
        
        
        # =========================================================
        # 【并行改造区】：同时启动独立的外部工具 (Step 1 & Step 4)
        # =========================================================
        print("\n🔍🔬 Step 1 & 4: Running MinCED and Prodigal concurrently...")
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            # 将两个耗时的外部命令同时丢进后台线程执行
            future_minced = executor.submit(self._detect_crispr_arrays, input_fasta, output_dir)
            future_prodigal = executor.submit(self._predict_proteins, input_fasta, output_dir)
            
            # 阻塞等待，直到两个任务都跑完，并获取它们的返回结果
            crispr_arrays = future_minced.result()
            protein_fasta = future_prodigal.result()
        # =========================================================

        # 检查 Step 1 (MinCED) 的结果
        if not crispr_arrays:
            print("⚠️  No CRISPR arrays detected!")
            return []
        
        # Step 2: 提取重复序列 (必须等 Step 1 跑完才能进行)
        print("\n🧬 Step 2: Extracting repeat sequences...")
        repeat_data = extract_repeat_sequences(input_fasta, crispr_arrays)
        if not repeat_data:
            print("⚠️  No repeat sequences extracted!")
            return []
        
        # Step 3: 重复序列分类
        print("\n🧠 Step 3: Classifying repeat sequences...")
        repeat_sequences = [r['sequence'] for r in repeat_data]
        print(f"repeat_sequences: {repeat_sequences}")
        repeat_predictions = self.repeat_classifier.predict(repeat_sequences)
        
        # 检查 Step 4 (Prodigal) 的结果 (实际上它已经跟 Step 1 一起跑完了)
        if not protein_fasta or not os.path.exists(protein_fasta):
            print("⚠️  Protein prediction failed!")
            protein_fasta = None
        
        cas_predictions = []
        hmm_results = []
        subtype_scores = {}
        
        if protein_fasta and os.path.exists(protein_fasta):
            # Step 5: Cas蛋白分类
            print("\n🧪 Step 5: Classifying Cas proteins...")
            protein_sequences = list(load_fasta_sequences(protein_fasta).values())
            protein_ids = list(load_fasta_sequences(protein_fasta).keys())
            if protein_sequences:
                cas_predictions = self.cas_classifier.predict(protein_sequences[:10])
                for i, pred in enumerate(cas_predictions):
                    if i < len(protein_ids):
                        pred['protein_id'] = protein_ids[i]
            
            # Step 6: 定向HMM扫描
            print("\n📊 Step 6: Directed HMM scanning with family-specific models...")
            # 【修复】：在这里提前初始化 protein_vectors，防止后面未定义
            protein_vectors = {} 
            if cas_predictions:
                hmm_results = self.hmm_scanner.scan_proteins_with_families(protein_fasta, cas_predictions)
                # 计算单蛋白倾向性向量
                protein_vectors = self.hmm_scanner.calculate_protein_propensity_vectors(hmm_results)
                print(f"✅ Calculated propensity vectors for {len(protein_vectors)} proteins")
        else:
            protein_vectors = {}

        # Step 7: 基于规则的最终分类
        print("\n🎯 Step 7: Final classification using rule-based approach...")
        final_results = self._integrate_results(
            repeat_data, repeat_predictions, cas_predictions, hmm_results, protein_vectors, output_dir
        )
        
        # 保存结果
        output_file = output_dir / f"{output_prefix}_results.csv"
        save_results_to_csv(final_results, str(output_file))
        
        print(f"\n✅ Pipeline completed successfully!")
        print(f"📄 Results saved to: {output_file}")
        
        return final_results
    
    def _detect_crispr_arrays(self, input_fasta: str, output_dir: Path) -> List[Dict]:
        """使用MinCED检测CRISPR阵列"""
        # MinCED 可以输出两个文件：常规输出和 GFF 格式输出。
        # 调用格式：java -jar minced.jar -gff <input.fa> <out.txt> <out.gff>
        minced_out = output_dir / "minced_output.out"
        minced_gff = output_dir / "minced_output.gff"
        minced_cmd = [
            "java", "-jar", MINCED_PATH,
            "-gff",
            input_fasta,
            str(minced_out),
            str(minced_gff)
        ]
        
        stdout, stderr, returncode = run_command(minced_cmd)
        
        if returncode != 0:
            print(f"MinCED failed: {stderr}")
            return []

        # 确保 GFF 文件已生成并包含 repeat_region
        if not os.path.exists(minced_gff):
            print(f"MinCED did not produce GFF output at {minced_gff}")
            return []

        return parse_minced_output(str(minced_gff))
    
    # def _predict_proteins(self, input_fasta: str, output_dir: Path) -> str:
    #     """使用Prodigal预测蛋白编码基因"""
    #     protein_fasta = output_dir / "predicted_proteins.faa"
    #     prodigal_cmd = [
    #         PRODIGAL_PATH,
    #         "-i", input_fasta,
    #         "-a", str(protein_fasta),
    #         "-p", "meta"  # 使用meta模式适用于宏基因组
    #     ]
        
    #     stdout, stderr, returncode = run_command(prodigal_cmd)
        
    #     if returncode != 0:
    #         print(f"Prodigal failed: {stderr}")
    #         return None
        
    #     return str(protein_fasta) if os.path.exists(protein_fasta) else None
    
    def _predict_proteins(self, input_fasta: str, output_dir: Path, threads: int = 32) -> str:
        """使用多进程切片加速版 Prodigal 预测蛋白编码基因"""
        import concurrent.futures
        import shutil
        
        protein_fasta = output_dir / "predicted_proteins.faa"
        chunk_dir = output_dir / "prodigal_chunks"
        chunk_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"      [Prodigal 加速] 正在将 FASTA 轮询分割为 {threads} 份以榨干多核算力...")
        
        chunk_paths = [chunk_dir / f"chunk_{i}.fna" for i in range(threads)]
        chunk_files = [open(path, 'w') for path in chunk_paths]
        
        try:
            with open(input_fasta, 'r') as f:
                idx = 0
                for line in f:
                    if line.startswith('>'):
                        idx = (idx + 1) % threads
                    chunk_files[idx].write(line)
        finally:
            for f in chunk_files:
                f.close()
                
        valid_chunks = [p for p in chunk_paths if os.path.getsize(p) > 0]
        
        # 【修改这里】：删掉了内部的 def run_prodigal_chunk
        
        print(f"      [Prodigal 加速] 启动 {len(valid_chunks)} 个并行任务轰炸...")
        faa_chunks = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            # 【修改这里】：调用刚才写在类外面的顶层函数 _run_prodigal_chunk_worker
            for faa_path in executor.map(_run_prodigal_chunk_worker, valid_chunks):
                if os.path.exists(faa_path) and os.path.getsize(faa_path) > 0:
                    faa_chunks.append(faa_path)
                    
        print("      [Prodigal 加速] 合并预测结果...")
        with open(protein_fasta, 'w') as outfile:
            for faa in faa_chunks:
                with open(faa, 'r') as infile:
                    outfile.write(infile.read())
                    
        shutil.rmtree(chunk_dir)
        
        return str(protein_fasta) if os.path.exists(protein_fasta) else None
    
    # 【修复】：彻底剔除旧的 subtype_scores 变量，改为接收 protein_vectors
    def _integrate_results(self, 
                          repeat_data: List[Dict],
                          repeat_predictions: List[Dict],
                          cas_predictions: List[Dict],
                          hmm_results: List[Dict],
                          protein_vectors: Dict[str, Dict[str, float]],
                          output_dir: Path) -> List[Dict]:
        """整合所有结果进行最终分类"""
        final_results = []
        cas_families = [pred['predicted_family'] for pred in cas_predictions]
        
        for i, repeat_info in enumerate(repeat_data):
            if i >= len(repeat_predictions):
                break
            repeat_pred = repeat_predictions[i]
            
            # 基于规则的最终分类 (聚合逻辑已在 rule_classifier 内部处理)
            rule_result = self.rule_classifier.classify_system(
                repeat_subtype=repeat_pred['predicted_subtype'],
                cas_families=cas_families,
                hmm_results=hmm_results,
                protein_vectors=protein_vectors  # 将向量传入决策层
            )
            
            result = {
                'array_id': repeat_info['array_id'],
                'sequence_id': repeat_info['seq_id'],
                'repeat_start': repeat_info['start'],
                'repeat_end': repeat_info['end'],
                'repeat_sequence': repeat_info['sequence'],
                'predicted_subtype': rule_result['predicted_subtype'],
                'confidence': rule_result['confidence'],
                'adaptation_score': rule_result['adaptation_score'],
                'interference_score': rule_result['interference_score'],
                'repeat_confidence': repeat_pred['confidence'],
                'detected_cas_families': ', '.join(rule_result['detected_cas_families']),
                'hmm_matches': ', '.join(rule_result['hmm_matches']) if rule_result['hmm_matches'] else 'None',
                'protein_vectors': str(protein_vectors) if protein_vectors else 'None'
            }
            
            final_results.append(result)
        
        return final_results