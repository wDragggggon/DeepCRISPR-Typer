import sys
import os

# 确保当前 `src` 目录在 sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from repeat_classifier import RepeatClassifier


if __name__ == '__main__':
    rc = RepeatClassifier()
    seqs = ['GTCGCGCCTTTACGGGCGCGTGGATTGAAAC', 'CGGTTCATCCCCACCTGCGTGGGGTTAAT']
    results = rc.predict(seqs)
    for r in results:
        print(r)
