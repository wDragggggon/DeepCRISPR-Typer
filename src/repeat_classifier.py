import torch
import torch.nn as nn
import numpy as np
from typing import List, Dict, Tuple
from pathlib import Path
import sys
import inspect
import os
import io
import importlib.util
import tempfile
import json
import subprocess
import runpy


cnn_att_source = "/home/wenlong/project/20260401_case9/工作代码/CRISPRclassify-CNN-Att/source"
if cnn_att_source not in sys.path:
    sys.path.append(cnn_att_source)


src_dir = str(Path(__file__).parent)
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)


project_root = str(Path(__file__).parent.parent)
if project_root not in sys.path:
    sys.path.insert(0, project_root)


try:
    from config.paths_config import CNN_ATT_MODEL_PATH
    from config.model_config import REPEAT_MODEL_CONFIG
except Exception:
    import importlib.util
    from pathlib import Path

    base_dir = Path(__file__).parent
    paths_config_path = base_dir / 'config' / 'paths_config.py'
    model_config_path = base_dir / 'config' / 'model_config.py'

    def load_module_from_path(name, path):
        spec = importlib.util.spec_from_file_location(name, str(path))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    if paths_config_path.exists():
        paths_mod = load_module_from_path('paths_config', paths_config_path)
        CNN_ATT_MODEL_PATH = getattr(paths_mod, 'CNN_ATT_MODEL_PATH')
    else:
        CNN_ATT_MODEL_PATH = None

    if model_config_path.exists():
        model_mod = load_module_from_path('model_config', model_config_path)
        REPEAT_MODEL_CONFIG = getattr(model_mod, 'REPEAT_MODEL_CONFIG')
    else:
        REPEAT_MODEL_CONFIG = {'num_subtypes': 33, 'max_repeat_length': 50}

base_dir = Path(__file__).parent
try:
    from src.utils import one_hot_encode_sequence
except Exception:
    try:
        from utils import one_hot_encode_sequence
    except Exception:
        utils_path = base_dir / 'utils.py'
        if utils_path.exists():
            utils_mod = load_module_from_path('utils', utils_path)
            one_hot_encode_sequence = getattr(utils_mod, 'one_hot_encode_sequence')
        else:

            def one_hot_encode_sequence(sequence: str, max_length: int = 50):
                nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 3}
                encoded = np.zeros((4, max_length), dtype=np.float32)
                sequence = sequence[:max_length].upper()
                for i, nucleotide in enumerate(sequence):
                    if nucleotide in nucleotide_to_index:
                        encoded[nucleotide_to_index[nucleotide], i] = 1.0
                return encoded

class RepeatClassifier:

    
    def __init__(self, model_path: str = None):
        self.model_path = model_path or CNN_ATT_MODEL_PATH
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model = self._load_model()
        self.subtype_names = self._get_subtype_names()
    
    def _load_model(self):
        # Diagnostic info
        print(f"RepeatClassifier: attempting to load model from: {self.model_path}")
        print(f"Type of model_path: {type(self.model_path)}")
        if isinstance(self.model_path, (str, bytes, os.PathLike)):
            print(f"Path exists: {os.path.exists(self.model_path)}")
        else:
            print("Model path is not a filesystem path")


        try:
            model_obj = torch.load(self.model_path, map_location=self.device)
            print("torch.load succeeded (direct)")
        except Exception as e1:
            print(f"torch.load(path) failed: {e1}")
            model_obj = None


        if model_obj is None and isinstance(self.model_path, (str, bytes, os.PathLike)) and os.path.exists(self.model_path):
            try:
                with open(self.model_path, 'rb') as f:
                    buf = io.BytesIO(f.read())
                    buf.seek(0)
                    model_obj = torch.load(buf, map_location=self.device)
                print("torch.load succeeded (from BytesIO)")
            except Exception as e2:
                print(f"torch.load(open(...)) failed: {e2}")
                model_obj = None


        if model_obj is None and isinstance(self.model_path, (str, bytes, os.PathLike)) and os.path.exists(self.model_path):
            try:
                model_obj = torch.jit.load(self.model_path, map_location=self.device)
                print("torch.jit.load succeeded")
            except Exception as e3:
                print(f"torch.jit.load failed: {e3}")
                model_obj = None


        if isinstance(model_obj, nn.Module):
            model_obj.eval()
            return model_obj.to(self.device)


        if isinstance(model_obj, dict):
            state_dict = None
            if 'model_state_dict' in model_obj:
                state_dict = model_obj['model_state_dict']
            elif 'state_dict' in model_obj:
                state_dict = model_obj['state_dict']
            else:

                state_dict = model_obj

            if state_dict is not None:

                try:
                    from CNN_Att import CNNClassifier

                    vocab_size = 5
                    embedding_dim = REPEAT_MODEL_CONFIG.get('embedding_dim', 64)
                    bio_feature_dim = getattr(self, 'bio_feature_dim', 2082)


                    inferred_num_classes = None
                    inferred_seq_length = None
                    conv_channels = REPEAT_MODEL_CONFIG.get('cnn_channels', 64)

                    for k, v in state_dict.items():
                        if k.endswith('fc_final.weight'):
                            try:
                                inferred_num_classes = int(v.shape[0])
                            except Exception:
                                pass
                        if k.endswith('fc_seq.weight'):
                            try:
                                in_features = int(v.shape[1])
                                if conv_channels and in_features % conv_channels == 0:
                                    inferred_seq_length = int(in_features // conv_channels + 9)
                            except Exception:
                                pass

                    num_classes = inferred_num_classes or REPEAT_MODEL_CONFIG.get('num_subtypes', 33)
                    seq_length = inferred_seq_length or REPEAT_MODEL_CONFIG.get('max_repeat_length', 48)

                    model = CNNClassifier(vocab_size, embedding_dim, int(num_classes), int(seq_length), bio_feature_dim)
                    try:

                        model.load_state_dict(state_dict)
                        model.eval()
                        print(f"Rebuilt CNNClassifier with seq_length={seq_length}, num_classes={num_classes} and loaded state_dict (strict=True)")
                        return model.to(self.device)
                    except Exception as e_load_strict:
                        print(f"Strict load failed: {e_load_strict} — trying non-strict load")
                        try:
                            model.load_state_dict(state_dict, strict=False)
                            model.eval()
                            print(f"Loaded state_dict with strict=False (seq_length={seq_length}, num_classes={num_classes})")
                            return model.to(self.device)
                        except Exception as e_load_loose:
                            print(f"Non-strict load also failed: {e_load_loose}")
                except Exception as e_imp:
                    print(f"Failed to import/rebuild CNNClassifier: {e_imp}")


        print("Using simplified simulation model for repeat classification")
        return self._create_simulation_model()

    def _call_crisprclassify(self, repeat_sequences: List[str]):
        """优先以模块方式调用 CRISPRclassify 的 infer 脚本；失败时以子进程调用。返回与 predict() 相同格式的结果列表或 None。"""

        from pathlib import Path

        if not repeat_sequences:
            return None

        tf = tempfile.NamedTemporaryFile(mode='w+', suffix='.fna', delete=False)
        fasta_path = Path(tf.name)
        try:
            for i, seq in enumerate(repeat_sequences):
                tf.write(f">repeat_{i}\n{seq}\n")
            tf.flush()
            tf.close()


            repo_infer = Path(__file__).resolve().parents[1] / 'crisprcastyper' / 'infer_fna.py'
            if not repo_infer.exists():
                return None

            try:
                spec = importlib.util.spec_from_file_location('dct_infer_fna', str(repo_infer))
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)


                if hasattr(mod, 'infer_from_fna') and callable(mod.infer_from_fna):
                    try:

                        result = mod.infer_from_fna(str(fasta_path), model_dir=None, device=self.device)
                    except TypeError:

                        result = mod.infer_from_fna(str(fasta_path), model_dir=None)


                    mapped = []
                    if isinstance(result, dict):
                        ids = result.get('ids') or []
                        seqs = result.get('sequences') or []
                        names = result.get('pred_names') or []
                        probs = result.get('pred_probs') or []
                        for i in range(len(seqs)):
                            mapped.append({
                                'sequence': seqs[i],
                                'predicted_subtype': names[i] if i < len(names) else None,
                                'confidence': float(max(probs[i]) if i < len(probs) else 0.0),
                                'probabilities': probs[i].tolist() if i < len(probs) else []
                            })
                        return mapped
                    elif isinstance(result, list):
                        return result
            except Exception as e:
                print(f"Failed to call crisprcastyper infer_fna: {e}")
                return None

            return None
        finally:
            try:
                fasta_path.unlink()
            except Exception:
                pass
    
    def _create_simulation_model(self):
        class SimpleRepeatClassifier(nn.Module):
            def __init__(self, num_subtypes=33):
                super().__init__()
                self.num_subtypes = num_subtypes
                self.feature_extractor = nn.Sequential(
                    nn.Conv1d(in_channels=4, out_channels=64, kernel_size=3, padding=1),
                    nn.ReLU(),
                    nn.AdaptiveAvgPool1d(1)
                )
                self.classifier = nn.Linear(64, num_subtypes)
            
            def forward(self, x):
                x = self.feature_extractor(x).squeeze(-1)
                logits = self.classifier(x)
                return torch.softmax(logits, dim=1)
        
        model = SimpleRepeatClassifier(REPEAT_MODEL_CONFIG['num_subtypes'])
        model.eval()
        return model.to(self.device)
    
    def _get_subtype_names(self):

        return [f"Type_{i+1}" for i in range(REPEAT_MODEL_CONFIG['num_subtypes'])]
    
    def predict(self, repeat_sequences: List[str]) -> List[Dict]:


        try:
            external_results = self._call_crisprclassify(repeat_sequences)
            if external_results:

                mapped = []
                if isinstance(external_results, dict):

                    for k, v in external_results.items():
                        mapped.append({
                            'sequence': k,
                            'predicted_subtype': v.get('predicted_subtype') or v.get('subtype') or v.get('label'),
                            'confidence': float(v.get('confidence', 0)),
                            'probabilities': v.get('probabilities', [])
                        })
                elif isinstance(external_results, list):
                    for item in external_results:
                        if isinstance(item, dict):
                            mapped.append({
                                'sequence': item.get('sequence') or item.get('seq') or '',
                                'predicted_subtype': item.get('predicted_subtype') or item.get('subtype') or item.get('label'),
                                'confidence': float(item.get('confidence', 0)),
                                'probabilities': item.get('probabilities', [])
                            })
                if mapped:
                    return mapped
        except Exception as e:
            print(f"External CRISPRclassify call failed: {e}")

        results = []

        for seq in repeat_sequences:

            encoded_seq = one_hot_encode_sequence(
                seq, 
                max_length=REPEAT_MODEL_CONFIG['max_repeat_length']
            )
            encoded_tensor = torch.tensor(encoded_seq, dtype=torch.float32).unsqueeze(0).to(self.device)
            

            with torch.no_grad():

                try:
                    sig = inspect.signature(self.model.forward)
                    params = list(sig.parameters.keys())
                except Exception:
                    params = []


                bio_dim = getattr(self.model, 'bio_feature_dim', REPEAT_MODEL_CONFIG.get('bio_feature_dim', 2082))
                bio_tensor = torch.zeros((encoded_tensor.size(0), bio_dim), dtype=torch.float32).to(self.device)
                try:
                    outputs = self.model(encoded_tensor, bio_tensor)
                except TypeError:
                    try:
                        outputs = self.model(encoded_tensor)
                    except Exception as e_call:
                        print(f"Model call failed with both signatures: {e_call}")
                        outputs = None


                if isinstance(outputs, torch.Tensor):
                    if outputs.dim() == 1:
                        outputs = outputs.unsqueeze(0)

                    probs = torch.softmax(outputs, dim=1).cpu().numpy()[0]
                else:
                    try:
                        probs = np.array(outputs)
                        if probs.ndim == 1:
                            probs = probs
                        else:
                            probs = probs[0]
                    except Exception:
                        probs = np.zeros(REPEAT_MODEL_CONFIG['num_subtypes'], dtype=float)

                probabilities = probs
                predicted_subtype_idx = int(np.argmax(probabilities))
                confidence = float(probabilities[predicted_subtype_idx])
                predicted_subtype = self.subtype_names[predicted_subtype_idx]
            
            results.append({
                'sequence': seq,
                'predicted_subtype': predicted_subtype,
                'confidence': confidence,
                'probabilities': probabilities.tolist()
            })
        
        return results