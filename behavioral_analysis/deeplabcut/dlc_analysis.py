import os
import glob
import shutil
import sys

import numpy as np
import pandas as pd
import cv2
import matplotlib

matplotlib.use("Agg")  # antes de importar pyplot/deeplabcut: não abre janelas, evita erro de thread
import matplotlib.pyplot as plt
import deeplabcut

# Para abrir a GUI (treino, rotulagem, avaliação):
#   conda activate DEEPLABCUT
#   QT_QPA_PLATFORM=xcb python -m deeplabcut

# =====================================================================
# CONFIGURAÇÃO
# =====================================================================
cfg         = '/home/arlindo/Desktop/tail-arlindo-2026-10-10/config.yaml'
project_dir = os.path.dirname(cfg)
videos_dir  = os.path.join(project_dir, 'videos')

EXTENSOES = ('.mp4', '.avi', '.mov', '.mkv', '.m4v')   # tipos aceitos (sem diferenciar maiúsculas)
EXCLUIR   = []   # nomes (sem extensão) a ignorar, ex.: ['MVI_0077', 'MVI_0078', 'MVI_0079']
SHUFFLE          = 2       # troque para 2 ao usar o modelo retreinado
FORCE_REANALYSIS = False   # True = reanalisa mesmo que o resultado já exista
JANELA_S = 300  # None = vídeo inteiro

# Parâmetros de desempenho (ajuste conforme a sua máquina)
BATCH_SIZE_ANALISE = 32      # suba para 64 se a VRAM permitir; volte um degrau se der out of memory
GERAR_VIDEO_ROTULADO = False  # False pula a etapa mais lenta quando só os dados interessam
DESENHAR_SKELETON = False    # menos desenho por frame no vídeo rotulado
TRAILPOINTS = 0              # 0 = sem rastro, mais rápido
GERAR_PLOTS = False           # False pula plot_trajectories

BODY = 'body'
NOSE = 'nose'
TAIL = 'tail'

PX_PER_CM    = 28.0   # estimativa grosseira; calibre com uma medida real
PCUTOFF      = 0.6    # descarta predições com likelihood menor
MAX_GAP_S    = 0.3    # interpola buracos de até 0,3 s
SMOOTH_S     = 0.2    # janela da mediana para reduzir tremor
THRESH_CM_S  = 2.0    # abaixo disso = parado (ponto de partida, precisa validar)
MIN_IMMOB_S  = 3.0    # duração mínima para contar como episódio de imobilidade

# ---------- pastas de saída ----------
work_dir = os.path.join(project_dir, '_trabalho')   # temporária (DeepLabCut escreve aqui)
pastas = {
    'tracking':            os.path.join(project_dir, 'tracking'),             # .h5, .csv, .pickle
    'videos_rotulados':    os.path.join(project_dir, 'videos_rotulados'),     # *_labeled.mp4
    'trajetorias':         os.path.join(project_dir, 'trajetorias'),          # plot_trajectories
    'graficos_velocidade': os.path.join(project_dir, 'graficos_velocidade'),  # PNGs de velocidade
    'dados':               os.path.join(project_dir, 'dados'),                # resumos em CSV
}
for p in [work_dir, *pastas.values()]:
    os.makedirs(p, exist_ok=True)


# =====================================================================
# LISTA DE VÍDEOS (vetor)
# =====================================================================
def listar_videos(pasta):
    caminhos = []
    for nome in sorted(os.listdir(pasta)):
        stem = os.path.splitext(nome)[0]
        caminho = os.path.join(pasta, nome)
        if not nome.lower().endswith(EXTENSOES):
            continue
        if 'DLC' in nome or stem in EXCLUIR:
            continue
        if not os.path.exists(caminho):          # link simbólico quebrado
            print(f'AVISO: ignorando {nome} (arquivo não encontrado, link quebrado?)')
            continue
        caminhos.append(caminho)
    return caminhos

videos = listar_videos(videos_dir)
if len(sys.argv) > 1:
    videos = [os.path.join(videos_dir, a) for a in sys.argv[1:]]

print(f'{len(videos)} vídeo(s) encontrado(s) em {videos_dir}:')

for v in videos:
    print('  ', os.path.basename(v))
if not videos:
    raise SystemExit('Nenhum vídeo para analisar.')

# =====================================================================
# FUNÇÕES
# =====================================================================
def achar_h5_filtrado(stem):
    """Retorna o .h5 filtrado mais recente deste vídeo e shuffle (ou None)."""
    padrao = os.path.join(pastas['tracking'], f'{stem}DLC*shuffle{SHUFFLE}*_filtered.h5')
    achados = glob.glob(padrao)
    return max(achados, key=os.path.getmtime) if achados else None

def rodar_deeplabcut(video, vtype, stem):
    """Analisa, filtra, gera vídeo rotulado e trajetórias; organiza os arquivos por tipo."""
    # 1. Análise: batch maior, o modelo é carregado uma vez por chamada
    deeplabcut.analyze_videos(
        cfg, [video], videotype=vtype, shuffle=SHUFFLE,
        save_as_csv=True, destfolder=work_dir,
        batch_size=BATCH_SIZE_ANALISE,
    )

    # 2. Filtro (rápido, roda sobre o .h5 já gerado)
    deeplabcut.filterpredictions(
        cfg, [video], videotype=vtype, shuffle=SHUFFLE,
        destfolder=work_dir,
    )

    # 3. Vídeo rotulado (a etapa mais lenta, feita na CPU)
    if GERAR_VIDEO_ROTULADO:
        deeplabcut.create_labeled_video(
            cfg, [video], videotype=vtype, shuffle=SHUFFLE,
            filtered=True, destfolder=work_dir,
            draw_skeleton=DESENHAR_SKELETON,
            trailpoints=TRAILPOINTS,
        )

    # 4. Trajetórias
    if GERAR_PLOTS:
        deeplabcut.plot_trajectories(
            cfg, [video], videotype=vtype, shuffle=SHUFFLE,
            filtered=True, destfolder=work_dir,
        )

    # 5. Organização dos arquivos (igual ao original)
    for f in glob.glob(os.path.join(work_dir, f'{stem}DLC*shuffle{SHUFFLE}*')):
        destino = pastas['videos_rotulados'] if f.endswith('_labeled.mp4') else pastas['tracking']
        shutil.move(f, os.path.join(destino, os.path.basename(f)))

    plots_origem = os.path.join(work_dir, 'plot-poses', stem)
    if os.path.isdir(plots_origem):
        plots_destino = os.path.join(pastas['trajetorias'], stem)
        shutil.rmtree(plots_destino, ignore_errors=True)  # sobrescreve a versão anterior
        shutil.move(plots_origem, plots_destino)

def speed_cm_s(df, fps, bp):
    ok = df[(bp, 'likelihood')] >= PCUTOFF
    x = df[(bp, 'x')].where(ok)
    y = df[(bp, 'y')].where(ok)
    gap = max(1, int(MAX_GAP_S * fps))
    x = x.interpolate(limit=gap, limit_area='inside')
    y = y.interpolate(limit=gap, limit_area='inside')
    win = max(1, int(SMOOTH_S * fps))
    x = x.rolling(win, center=True, min_periods=1).median()
    y = y.rolling(win, center=True, min_periods=1).median()
    step_px = np.hypot(x.diff(), y.diff())
    return step_px / PX_PER_CM * fps, ok.mean()  # cm/s e fração de frames confiáveis


def bouts(mask, min_s, fps):
    m = np.asarray(mask, dtype=int)
    d = np.diff(np.r_[0, m, 0])
    starts, ends = np.where(d == 1)[0], np.where(d == -1)[0]
    return [(s, e) for s, e in zip(starts, ends) if (e - s) / fps >= min_s]


def calcular_metricas(h5, video, stem):
    df = pd.read_hdf(h5)
    df.columns = df.columns.droplevel(0)  # remove o nível "scorer"
    bodyparts = list(df.columns.get_level_values(0).unique())

    faltando = [bp for bp in (BODY, NOSE, TAIL) if bp not in bodyparts]
    assert not faltando, f'Pontos não encontrados: {faltando}. Disponíveis: {bodyparts}'

    cap = cv2.VideoCapture(video)
    fps = cap.get(cv2.CAP_PROP_FPS)
    largura = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    altura = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    cap.release()
    assert fps and fps > 0, f'Não foi possível ler o fps de {video}'

    if JANELA_S:
        df = df.iloc[:int(round(JANELA_S * fps))]   # 300 s * 29,97 fps = 8991 frames

    total = len(df) / fps
    print(f'  fps = {fps:.2f}, frames = {len(df)}, duração = {total:.1f} s')

    rows = {}
    for name, bp in [('corpo', BODY), ('nariz', NOSE), ('cauda', TAIL)]:
        spd, frac_ok = speed_cm_s(df, fps, bp)
        immobile = (spd < THRESH_CM_S).fillna(False)
        ep = bouts(immobile, MIN_IMMOB_S, fps)
        t_immob = sum(e - s for s, e in ep) / fps
        rows[name] = {
            'frames_confiaveis_%': round(100 * frac_ok, 1),
            'distancia_total_cm': round((spd / fps).sum(), 1),
            'velocidade_media_cm_s': round(spd.mean(), 2),
            'tempo_parado_s': round(t_immob, 1),
            'tempo_movimento_s': round(total - t_immob, 1),
            'n_episodios_imobilidade': len(ep),
            'latencia_1a_imobilidade_s': round(ep[0][0] / fps, 1) if ep else None,
        }

        t = np.arange(len(spd)) / fps
        plt.figure(figsize=(12, 3))
        plt.plot(t, spd, lw=0.6)
        for s, e in ep:
            plt.axvspan(s / fps, e / fps, color='orange', alpha=0.3)
        plt.axhline(THRESH_CM_S, color='r', ls='--', lw=0.8)
        plt.xlabel('tempo (s)')
        plt.ylabel('velocidade (cm/s)')
        plt.title(f'{stem} - {name}')
        plt.tight_layout()
        plt.savefig(os.path.join(pastas['graficos_velocidade'],
                                 f'{stem}_{name}_velocidade_shuffle{SHUFFLE}.png'), dpi=150)
        plt.close()

    # imobilidade simultânea de TODOS os pontos
    speeds = {bp: speed_cm_s(df, fps, bp)[0] for bp in bodyparts}
    imm_por_ponto = pd.DataFrame({bp: (s < THRESH_CM_S) for bp, s in speeds.items()})
    todos_parados = imm_por_ponto.all(axis=1)  # True só se todos os pontos estão parados
    ep_all = bouts(todos_parados, MIN_IMMOB_S, fps)
    t_all = sum(e - s for s, e in ep_all) / fps
    validos_all = pd.DataFrame(speeds).notna().all(axis=1).mean()

    rows['todos_pontos'] = {
        'frames_confiaveis_%': round(100 * validos_all, 1),
        'tempo_parado_s': round(t_all, 1),
        'tempo_movimento_s': round(total - t_all, 1),
        'n_episodios_imobilidade': len(ep_all),
        'latencia_1a_imobilidade_s': round(ep_all[0][0] / fps, 1) if ep_all else None,
    }
    resumo = pd.DataFrame(rows).T.infer_objects()

    # metadados (ficam em todas as linhas)
    cap = cv2.VideoCapture(video)
    largura = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    altura = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    cap.release()
    treino = set(os.listdir(os.path.join(project_dir, 'labeled-data')))

    resumo.insert(0, 'no_treino', stem in treino)
    resumo.insert(1, 'resolucao', f'{largura}x{altura}')
    resumo.insert(2, 'fps', round(fps, 2))
    resumo.insert(3, 'janela_s', round(total, 1))
    resumo.insert(4, 'shuffle', SHUFFLE)
    resumo.insert(5, 'limiar_cm_s', THRESH_CM_S)
    resumo.insert(6, 'min_imob_s', MIN_IMMOB_S)
    resumo.insert(7, 'px_por_cm', PX_PER_CM)
    resumo['baixa_confianca'] = resumo['frames_confiaveis_%'] < 90  # revisar esses vídeos
    return resumo


# =====================================================================
# LOOP SOBRE OS VÍDEOS
# =====================================================================
resultados = {}  # stem -> tabela de métricas
falhas = {}  # stem -> mensagem de erro

for i, video in enumerate(videos, start=1):
    nome = os.path.basename(video)
    stem = os.path.splitext(nome)[0]
    vtype = os.path.splitext(nome)[1]
    print(f'\n[{i}/{len(videos)}] {nome}')

    try:
        h5 = achar_h5_filtrado(stem)
        if h5 is None or FORCE_REANALYSIS:
            rodar_deeplabcut(video, vtype, stem)
            h5 = achar_h5_filtrado(stem)
        else:
            print(f'  Resultado já existe, pulando análise: {os.path.basename(h5)}')

        if h5 is None:
            raise FileNotFoundError(
                f'Nenhum .h5 filtrado para {stem} (shuffle {SHUFFLE}) em {pastas["tracking"]}')

        resumo = calcular_metricas(h5, video, stem)
        print(resumo)
        resumo.to_csv(os.path.join(pastas['dados'],
                                   f'{stem}_resumo_movimento_shuffle{SHUFFLE}.csv'))
        resultados[stem] = resumo

    except Exception as e:                      # um vídeo com problema não derruba o lote
        falhas[stem] = f'{type(e).__name__}: {e}'
        print(f'  FALHOU: {falhas[stem]}')

# =====================================================================
# RESUMO GERAL (lido dos CSVs individuais em dados/)
# =====================================================================
partes = []
for f in sorted(glob.glob(os.path.join(pastas['dados'], f'MVI_*_resumo_movimento_shuffle{SHUFFLE}.csv'))):
    t = pd.read_csv(f, index_col=0)
    t.index.name = 'ponto'
    t.insert(0, 'video', os.path.basename(f).split('_resumo_movimento')[0])
    partes.append(t.reset_index())

if partes:
    geral = pd.concat(partes, ignore_index=True)
    caminho_geral = os.path.join(pastas['dados'], f'resumo_geral_shuffle{SHUFFLE}.csv')
    geral.to_csv(caminho_geral, index=False)
    print(f'\nResumo geral ({len(partes)} vídeos) salvo em: {caminho_geral}')

print(f'\nConcluído: {len(resultados)} ok, {len(falhas)} com falha.')
for stem, msg in falhas.items():
    print(f'  - {stem}: {msg}')

