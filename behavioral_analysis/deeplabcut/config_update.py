import os
import shutil
import cv2
from deeplabcut.utils import auxiliaryfunctions

cfg        = '/home/arlindo/Desktop/tail-arlindo-2026-10-10/config.yaml'
treino_dir = os.path.join(os.path.dirname(cfg), 'treino')
EXTENSOES  = ('.mp4', '.avi', '.mov', '.mkv', '.m4v')

# backup do config antes de mexer
shutil.copy(cfg, cfg + '.bak')

cfg_dict = auxiliaryfunctions.read_config(cfg)

video_sets = {}
for nome in sorted(os.listdir(treino_dir)):
    if not nome.lower().endswith(EXTENSOES):
        continue
    caminho = os.path.join(treino_dir, nome)
    cap = cv2.VideoCapture(caminho)
    largura = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    altura  = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    cap.release()
    if largura == 0 or altura == 0:
        print(f'AVISO: não consegui ler {nome}, ignorando')
        continue
    video_sets[caminho] = {'crop': f'0, {largura}, 0, {altura}'}

cfg_dict['video_sets'] = video_sets
auxiliaryfunctions.write_config(cfg, cfg_dict)

print(f'{len(video_sets)} vídeo(s) em video_sets:')
for v in video_sets:
    print('  ', v)