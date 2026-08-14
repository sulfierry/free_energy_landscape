import os
import sys
import shutil
import hashlib
import platform
import tempfile
import subprocess
import numpy as np
import multiprocessing
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import proj3d
from scipy import ndimage
from scipy.stats import gaussian_kde
from joblib import Parallel, delayed
from matplotlib.colors import LinearSegmentedColormap

class FreeEnergyLandscape:

    def __init__(self, cv1_path, cv2_path,
                 temperature, boltzmann_constant,
                 bins=100, kde_bandwidth=None,
                 cv_names=['CV1', 'CV2'], discrete=None,
                 xlim_inf=None, xlim_sup=None, ylim_inf=None, ylim_sup=None,
                 max_energy=None, grid_size=100, n_jobs=None,
                 basin_method='watershed', basin_connectivity=2,
                 basin_min_frames=1, basin_min_depth=None,
                 mask_min_count=1.0, cmap='viridis_r', dpi=200):


        self.cv1_path = cv1_path
        self.cv2_path = cv2_path
        self.temperature = temperature
        self.kB = boltzmann_constant
        self.cv_names = cv_names
        self.colors = [
            (0.0, "darkblue"),
            (0.1, "blue"),
            (0.2, "dodgerblue"),
            (0.3, "deepskyblue"),
            (0.4, "lightblue"),
            (0.5, "azure"),
            (0.6, "oldlace"),
            (0.7, "wheat"),
            (0.8, "lightcoral"),
            (0.9, "indianred"),
            (1.0, "darkred")
        ]
        self.custom_cmap = LinearSegmentedColormap.from_list(
            "custom_energy",
            self.colors
            )


        self.proj1_data_original = None
        self.proj2_data_original = None
        self.frames = None
        self.bins = bins
        self.kde_bandwidth = kde_bandwidth
        self.discrete = discrete
        self.max_energy = max_energy
        self.grid_size = int(grid_size)
        self.n_jobs = n_jobs if n_jobs else multiprocessing.cpu_count()
        self.basin_method = basin_method
        self.basin_connectivity = int(basin_connectivity)
        self.basin_min_frames = int(basin_min_frames)
        # Profundidade minima para uma bacia contar como estado proprio. 'auto' escolhe
        # o corte no maior degrau do espectro de persistencia, limitado a kB*T: acima da
        # escala termica a barreira e sempre real, abaixo do degrau e ruido do KDE.
        self.basin_min_depth = ('auto' if basin_min_depth is None
                                else float(basin_min_depth))
        # Uma celula da grade so conta como amostrada se o KDE preve pelo menos
        # esse numero de frames dentro dela. Evita extrapolar para o vazio.
        self.mask_min_count = float(mask_min_count)
        self.cmap_name = cmap
        self.dpi = int(dpi)
        self._basins = None

        # Kernel unico, reutilizado por todos os caminhos de calculo. A chave garante
        # que o cache seja invalidado se os dados ou a largura de banda mudarem.
        self._kernel = None
        self._kernel_key = None
        self._cache_key = None
        self.cached_results = None
        self._frame_energy = None

        self.discreet_colors = [
            'purple', 
            'darkorange', 
            'green', 
            'lightgrey', 
            'red',
            'magenta',
            'mediumorchid',
            'deeppink',
            'peru',
            'indianred'
            ]
    
        self.discreet_markers = [
            '*', 
            's', 
            '^', 
            'D', 
            'o',
            'p',
            'h',
            'v',
            'X',
            'd'
            ]

        # Novos atributos para limites dos eixos
        self.xlim_inf = xlim_inf
        self.xlim_sup = xlim_sup
        self.ylim_inf = ylim_inf
        self.ylim_sup = ylim_sup

    def load_data(self):
        """Carrega as CVs e os indices de frame, validando a consistencia entre os arquivos."""
        raw1 = np.loadtxt(self.cv1_path)
        raw2 = np.loadtxt(self.cv2_path)

        for path, raw in ((self.cv1_path, raw1), (self.cv2_path, raw2)):
            if raw.ndim != 2 or raw.shape[1] < 2:
                raise ValueError(
                    f"'{path}' must have at least 2 columns (frame index and CV value); "
                    f"got shape {raw.shape}."
                )

        if raw1.shape[0] != raw2.shape[0]:
            raise ValueError(
                f"CV files have different numbers of frames: "
                f"'{self.cv1_path}' has {raw1.shape[0]}, '{self.cv2_path}' has {raw2.shape[0]}."
            )

        frames1 = raw1[:, 0].astype(np.int64)
        frames2 = raw2[:, 0].astype(np.int64)
        if not np.array_equal(frames1, frames2):
            n_bad = int(np.count_nonzero(frames1 != frames2))
            first_bad = int(np.flatnonzero(frames1 != frames2)[0])
            msg = [
                f"Frame indices do not match between '{self.cv1_path}' and '{self.cv2_path}' "
                f"({n_bad} mismatching rows, first at row {first_bad}). "
                f"Row {first_bad} is frame {frames1[first_bad]} in CV1 but frame "
                f"{frames2[first_bad]} in CV2, so the two CVs would be paired across "
                f"different conformations."
            ]
            # Diagnostico comum: trajetorias concatenadas cujo contador reinicia em
            # linhas diferentes nos dois arquivos.
            r1 = np.flatnonzero(np.diff(frames1) < 0) + 1
            r2 = np.flatnonzero(np.diff(frames2) < 0) + 1
            if (len(r1) or len(r2)) and not np.array_equal(r1, r2):
                msg.append(
                    f"The frame counter restarts at rows {r1.tolist()} in CV1 and "
                    f"{r2.tolist()} in CV2: the files look like concatenated trajectories "
                    f"whose blocks have different lengths. Trim the blocks so both files "
                    f"restart at the same rows before running the analysis."
                )
            raise ValueError(" ".join(msg))

        self.frames = frames1
        self.proj1_data_original = raw1[:, 1]
        self.proj2_data_original = raw2[:, 1]

        for name, values in ((self.cv_names[0], self.proj1_data_original),
                             (self.cv_names[1], self.proj2_data_original)):
            if not np.all(np.isfinite(values)):
                raise ValueError(f"CV '{name}' contains NaN or infinite values.")

        # Dados novos invalidam qualquer kernel/superficie em cache
        self._kernel = None
        self._kernel_key = None
        self._cache_key = None
        self.cached_results = None
        self._frame_energy = None

    @staticmethod
    def _array_key(data):
        """Chave estavel para invalidacao de cache."""
        buf = np.ascontiguousarray(data, dtype=np.float64)
        digest = hashlib.blake2b(buf.tobytes(), digest_size=16).hexdigest()
        return (buf.shape, digest)

    def get_kernel(self, data):
        """Ajusta (uma unica vez) um gaussian_kde sobre TODOS os dados.

        Um kernel unico e a unica forma correta de comparar densidades entre pontos.
        Ajustar um kernel por bloco da trajetoria torna as densidades incomparaveis
        e faz o resultado depender do numero de CPUs.
        """
        key = (self._array_key(data), self.kde_bandwidth)
        if self._kernel is None or self._kernel_key != key:
            self._kernel = gaussian_kde(data.T, bw_method=self.kde_bandwidth)
            self._kernel_key = key
        return self._kernel

    def evaluate_density(self, kernel, points, n_jobs=None):
        """Avalia um kernel JA ajustado sobre `points` (N, 2), paralelizando so a avaliacao.

        O particionamento afeta apenas a divisao do trabalho: a concatenacao dos
        chunks e identica, ponto a ponto, a `kernel(points.T)` em serie.
        """
        points = np.ascontiguousarray(points, dtype=np.float64)
        n_jobs = n_jobs if n_jobs else self.n_jobs
        n_chunks = max(1, min(int(n_jobs), points.shape[0]))
        if n_chunks == 1:
            return kernel(points.T)
        chunks = np.array_split(points, n_chunks, axis=0)
        densities = Parallel(n_jobs=n_chunks)(
            delayed(kernel)(chunk.T) for chunk in chunks
        )
        return np.concatenate(densities)

    def _to_free_energy(self, density, cap=None, mask=None):
        """Converte densidade em energia livre relativa: G = -kB*T*ln(rho), com min em 0.

        Onde `mask` e False a celula nao foi amostrada e recebe NaN: o KDE ali esta
        extrapolando para regiao sem um unico frame, e pintar esse vazio como se fosse
        energia medida domina a figura e distorce a escala de cor.
        """
        density = np.clip(np.asarray(density, dtype=np.float64),
                          np.finfo(np.float64).tiny, None)
        G = -self.kB * self.temperature * np.log(density)
        if mask is not None:
            G = np.where(mask, G, np.nan)
        G = G - np.nanmin(G)
        if cap is not None:
            G = np.clip(G, 0, cap)
        return G

    def get_cmap(self):
        """Colormap da superficie. 'classic' devolve o mapa azul-vermelho das versoes 1.x."""
        if self.cmap_name == 'classic':
            return self.custom_cmap
        return plt.get_cmap(self.cmap_name)

    def calculate_free_energy(self, data):
        """Superficie de energia livre numa grade regular, a partir do kernel global."""
        data = np.ascontiguousarray(data, dtype=np.float64)

        # Ajusta a geração da grade para respeitar os limites especificados
        x_min = self.xlim_inf if self.xlim_inf is not None else data[:, 0].min()
        x_max = self.xlim_sup if self.xlim_sup is not None else data[:, 0].max()
        y_min = self.ylim_inf if self.ylim_inf is not None else data[:, 1].min()
        y_max = self.ylim_sup if self.ylim_sup is not None else data[:, 1].max()

        key = (self._array_key(data), self.kde_bandwidth, self.grid_size,
               self.max_energy, self.mask_min_count, x_min, x_max, y_min, y_max)
        if self.cached_results is not None and self._cache_key == key:
            return self.cached_results

        kernel = self.get_kernel(data)

        n = self.grid_size
        X_original, Y_original = np.mgrid[x_min:x_max:n * 1j, y_min:y_max:n * 1j]
        positions_original = np.column_stack([X_original.ravel(), Y_original.ravel()])
        Z_original = np.reshape(self.evaluate_density(kernel, positions_original),
                                X_original.shape)

        # Numero esperado de frames dentro de cada celula da grade. Celulas onde nem
        # um frame e esperado nao foram amostradas: o KDE ali e extrapolacao pura.
        cell_area = ((x_max - x_min) / (n - 1)) * ((y_max - y_min) / (n - 1))
        expected_counts = Z_original * cell_area * data.shape[0]
        mask = expected_counts >= self.mask_min_count
        # A borda bruta fica rendilhada e fragmenta a superficie 3D em dentes de serra.
        # Fechamento morfologico + preenchimento de buracos deixam um contorno continuo
        # sem mexer no interior da regiao amostrada.
        mask = ndimage.binary_fill_holes(
            ndimage.binary_closing(mask, structure=np.ones((3, 3), dtype=bool)))
        if not mask.any():
            raise ValueError(
                f"No grid cell reaches the sampling threshold of {self.mask_min_count} "
                f"expected frames. Lower --mask_min_count or coarsen --grid."
            )

        G_original = self._to_free_energy(Z_original, cap=self.max_energy, mask=mask)
        # Mesma superficie sem o recorte, usada so para desenhar o 3D sem buracos
        G_unmasked = self._to_free_energy(Z_original, cap=None)
        G_unmasked = G_unmasked - (np.nanmin(np.where(mask, G_unmasked, np.nan)))

        self.cached_results = {'X_original': X_original,
                            'Y_original': Y_original,
                            'G_original': G_original,
                            'G_unmasked': G_unmasked,
                            'mask': mask,
                            'expected_counts': expected_counts
                            }
        self._cache_key = key
        return self.cached_results


    def combined_data(self):
        """Matriz (N, 2) com as duas CVs por frame."""
        if self.proj1_data_original is None or self.proj2_data_original is None:
            raise ValueError("Data not loaded. Run load_data first.")
        return np.column_stack([self.proj1_data_original, self.proj2_data_original])

    def frame_energies(self):
        """Energia livre de cada frame amostrado (kJ/mol, minimo em 0). Calculada uma vez."""
        if self._frame_energy is None:
            data = self.combined_data()
            kernel = self.get_kernel(data)
            self._frame_energy = self._to_free_energy(self.evaluate_density(kernel, data))
        return self._frame_energy

    def _grid_indices(self, result):
        """Celula da grade que contem cada frame. Retorna (i, j, dentro_dos_limites)."""
        X, Y = result['X_original'], result['Y_original']
        x_axis, y_axis = X[:, 0], Y[0, :]
        nx, ny = x_axis.size, y_axis.size
        dx = (x_axis[-1] - x_axis[0]) / (nx - 1)
        dy = (y_axis[-1] - y_axis[0]) / (ny - 1)

        i = np.rint((self.proj1_data_original - x_axis[0]) / dx).astype(np.int64)
        j = np.rint((self.proj2_data_original - y_axis[0]) / dy).astype(np.int64)
        inside = (i >= 0) & (i < nx) & (j >= 0) & (j < ny)
        return np.clip(i, 0, nx - 1), np.clip(j, 0, ny - 1), inside

    @staticmethod
    def _steepest_descent_labels(G):
        """Watershed por descida de gradiente na grade.

        Cada celula aponta para o vizinho (entre os 8) de menor energia. Seguindo esses
        ponteiros ate o ponto fixo, toda celula chega a um minimo local: a sua bacia de
        atracao. Celulas NaN (fora da regiao amostrada) ficam com rotulo 0.

        Nao existe ciclo de tamanho 2 porque o passo so aceita vizinho estritamente
        menor, entao a iteracao de ponteiros converge sempre.
        """
        n_rows, n_cols = G.shape
        valid = np.isfinite(G)
        walls = np.where(valid, G, np.inf)

        flat = np.arange(G.size).reshape(G.shape)
        target = flat.copy()
        best = walls.copy()

        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                shifted = np.full_like(walls, np.inf)
                idx = np.zeros_like(flat)
                src_i = slice(max(0, -di), n_rows - max(0, di))
                dst_i = slice(max(0, di), n_rows - max(0, -di))
                src_j = slice(max(0, -dj), n_cols - max(0, dj))
                dst_j = slice(max(0, dj), n_cols - max(0, -dj))
                shifted[src_i, src_j] = walls[dst_i, dst_j]
                idx[src_i, src_j] = flat[dst_i, dst_j]

                downhill = shifted < best
                best = np.where(downhill, shifted, best)
                target = np.where(downhill, idx, target)

        # Iteracao de ponteiros ate o ponto fixo (minimos apontam para si mesmos)
        parent = target.ravel().copy()
        while True:
            jumped = parent[parent]
            if np.array_equal(jumped, parent):
                break
            parent = jumped

        roots = np.unique(parent[valid.ravel()])
        # Bacia 1 e a mais profunda
        roots = roots[np.argsort(G.ravel()[roots])]
        lookup = np.zeros(G.size, dtype=np.int64)
        lookup[roots] = np.arange(1, roots.size + 1)

        labels = lookup[parent].reshape(G.shape)
        labels[~valid] = 0
        return labels

    @staticmethod
    def _basin_saddles(G, labels):
        """Menor sela entre cada par de bacias adjacentes.

        A sela entre A e B e o menor `max(G[x], G[y])` sobre todos os pares de celulas
        vizinhas com x em A e y em B: o ponto mais baixo da crista que as separa.
        """
        saddles = {}
        # 4 deslocamentos cobrem toda a vizinhanca de 8 sem repetir par
        for di, dj in ((0, 1), (1, 0), (1, 1), (1, -1)):
            a = labels[max(0, -di):labels.shape[0] - max(0, di),
                       max(0, -dj):labels.shape[1] - max(0, dj)]
            b = labels[max(0, di):labels.shape[0] - max(0, -di),
                       max(0, dj):labels.shape[1] - max(0, -dj)]
            ga = G[max(0, -di):G.shape[0] - max(0, di),
                   max(0, -dj):G.shape[1] - max(0, dj)]
            gb = G[max(0, di):G.shape[0] - max(0, -di),
                   max(0, dj):G.shape[1] - max(0, -dj)]

            border = (a != b) & (a > 0) & (b > 0)
            if not border.any():
                continue
            heights = np.maximum(ga[border], gb[border])
            for la, lb, h in zip(a[border], b[border], heights):
                pair = (int(min(la, lb)), int(max(la, lb)))
                if h < saddles.get(pair, np.inf):
                    saddles[pair] = float(h)
        return saddles

    def _auto_min_depth(self, persistence):
        """Corte de profundidade tirado do maior degrau do espectro de persistencia.

        Ordenadas as profundidades das bacias, estrutura real e ruido do KDE costumam
        ficar separados por um salto largo. Cortar nesse salto e mais robusto que fixar
        um valor: uma escolha rigida ou funde estados verdadeiros (quando todas as
        barreiras do sistema sao rasas) ou deixa passar rugosidade do estimador.

        O corte nunca passa de kB*T: acima da escala termica a barreira e real por
        definicao e nao ha motivo para fundir.
        """
        thermal = self.kB * self.temperature
        if len(persistence) < 2:
            return 0.0
        gaps = [(persistence[i] - persistence[i + 1], i) for i in range(len(persistence) - 1)]
        gap, i = max(gaps)
        if gap <= 0:
            return 0.0
        cut = 0.5 * (persistence[i] + persistence[i + 1])
        return float(min(cut, thermal))

    def _merge_shallow_basins(self, G, labels, min_depth):
        """Funde bacias cuja profundidade (prominencia topologica) fica abaixo de min_depth.

        Profundidade de uma bacia = altura da sela mais baixa que a separa de uma vizinha,
        menos a energia do seu proprio minimo. Uma barreira menor que kB*T e cruzada por
        agitacao termica, entao os dois lados nao sao estados distintos: o KDE apenas
        rugou a superficie. Fundir por profundidade e o criterio topologico padrao
        (persistencia) e nao depende da resolucao da grade.
        """
        labels = labels.copy()
        if min_depth <= 0:
            return labels

        while True:
            active = np.unique(labels[labels > 0])
            if active.size <= 1:
                break

            minima = {int(lab): float(np.nanmin(np.where(labels == lab, G, np.nan)))
                      for lab in active}
            saddles = self._basin_saddles(G, labels)
            if not saddles:
                break

            shallowest, victim, winner = np.inf, None, None
            for (a, b), height in saddles.items():
                for this, other in ((a, b), (b, a)):
                    depth = height - minima[this]
                    # A bacia mais rasa do par e a candidata a ser absorvida
                    if minima[this] >= minima[other] and depth < shallowest:
                        shallowest, victim, winner = depth, this, other

            if victim is None or shallowest >= min_depth:
                break
            labels[labels == victim] = winner

        # Renumera 1..N por energia crescente do minimo
        active = np.unique(labels[labels > 0])
        order = sorted(active, key=lambda lab: np.nanmin(np.where(labels == lab, G, np.nan)))
        lookup = np.zeros(int(labels.max()) + 1, dtype=np.int64)
        for new_id, lab in enumerate(order, start=1):
            lookup[lab] = new_id
        return lookup[labels]

    def identify_basins(self, threshold=None, method=None, connectivity=2,
                        min_frames=1, min_depth=None):
        """Agrupa a superficie em bacias e devolve um ponto representativo de cada.

        method='watershed' (padrao)
            Bacia de atracao: cada celula da grade desce pelo gradiente ate um minimo
            local, e o minimo alcancado define a bacia. Particiona toda a regiao
            amostrada, nao precisa de threshold e nao deixa frame orfao. Bacias
            separadas por barreiras mais rasas que `min_depth` sao fundidas.

        method='connected'
            Componentes conexos de {G <= threshold}. Simples, mas so agrupa o que esta
            abaixo do corte e a resposta muda com o corte escolhido.

        `threshold`, no modo watershed, e apenas um filtro de exibicao: so as bacias cujo
        minimo fica abaixo dele sao reportadas.
        """
        method = method or self.basin_method
        if method == 'connected':
            if threshold is None:
                raise ValueError("method='connected' requires an energy threshold (--energy).")
            return self._identify_basins_connected(threshold, connectivity, min_frames)

        if method != 'watershed':
            raise ValueError(f"Unknown basin method {method!r}; use 'watershed' or 'connected'.")

        min_depth = self.basin_min_depth if min_depth is None else float(min_depth)

        data = self.combined_data()

        result = self.calculate_free_energy(data)
        G, X, Y = result['G_original'], result['X_original'], result['Y_original']

        raw_labels = self._steepest_descent_labels(G)
        # Espectro de persistencia antes da fusao: mostra o que min_depth esta engolindo
        raw_saddles = self._basin_saddles(G, raw_labels)
        raw_minima = {int(lab): float(np.nanmin(np.where(raw_labels == lab, G, np.nan)))
                      for lab in np.unique(raw_labels[raw_labels > 0])}
        persistence = sorted(
            (min(h for (a, b), h in raw_saddles.items() if lab in (a, b)) - g
             for lab, g in raw_minima.items()
             if any(lab in pair for pair in raw_saddles)),
            reverse=True)

        auto_depth = min_depth == 'auto'
        if auto_depth:
            min_depth = self._auto_min_depth(persistence)

        labels = self._merge_shallow_basins(G, raw_labels, min_depth)

        energies = self.frame_energies()
        i, j, inside = self._grid_indices(result)
        frame_label = np.where(inside, labels[i, j], 0)
        saddles = self._basin_saddles(G, labels)

        candidates = []
        for lab in np.unique(labels[labels > 0]):
            lab = int(lab)
            cells = labels == lab
            members = np.flatnonzero(frame_label == lab)
            if members.size < min_frames:
                continue
            gi, gj = np.unravel_index(np.nanargmin(np.where(cells, G, np.nan)), G.shape)
            g_min = float(G[gi, gj])
            if threshold is not None and g_min > threshold:
                continue
            neighbours = [h for (a, b), h in saddles.items() if lab in (a, b)]
            rep = members[np.argmin(energies[members])]
            candidates.append({
                'G_min': g_min,
                'depth': float(min(neighbours) - g_min) if neighbours else float('inf'),
                'cv1_min': float(X[gi, gj]),
                'cv2_min': float(Y[gi, gj]),
                'n_frames': int(members.size),
                'population': float(members.size / len(energies)),
                'rep_frame': int(self.frames[rep]),
                'rep_cv1': float(self.proj1_data_original[rep]),
                'rep_cv2': float(self.proj2_data_original[rep]),
                'rep_energy': float(energies[rep]),
                'area_cells': int(cells.sum()),
                '_label': lab,
            })

        candidates.sort(key=lambda b: b['G_min'])

        shown = np.zeros_like(labels)
        frame_basin = np.zeros_like(frame_label)
        for new_id, basin in enumerate(candidates, start=1):
            shown[labels == basin['_label']] = new_id
            frame_basin[frame_label == basin['_label']] = new_id
            basin['id'] = new_id
            del basin['_label']

        return {'basins': candidates, 'labels': shown, 'frame_basin': frame_basin,
                'threshold': threshold, 'method': 'watershed',
                'persistence': persistence, 'min_depth': min_depth,
                'auto_depth': auto_depth, 'n_raw_minima': len(raw_minima)}

    def _identify_basins_connected(self, threshold, connectivity=2, min_frames=1):
        """Agrupa a regiao de baixa energia em bacias e devolve um ponto representativo de cada.

        Uma bacia e um componente conexo de {G <= threshold} na superficie de energia livre.
        Duas regioes separadas por uma barreira acima do threshold sao bacias distintas, que e
        a definicao operacional de estado metaestavel nessa escala de energia. Isso substitui
        a marcacao ponto-a-ponto por faixa de energia, que coloria nos da grade e nao
        conformacoes: dois pontos com a mesma energia podiam estar em bacias sem conexao.

        Parametros
        ----------
        threshold : float
            Corte de energia (kJ/mol) que define a regiao agrupada.
        connectivity : int
            1 = vizinhanca de 4 celulas, 2 = vizinhanca de 8 (padrao, menos sensivel a
            artefatos de discretizacao da grade).
        min_frames : int
            Bacias com menos frames amostrados que isso sao descartadas como artefato do KDE.

        Retorna
        -------
        dict com 'basins' (lista ordenada por energia crescente), 'labels' (mapa da grade
        com o id da bacia por celula, 0 = fora) e 'frame_basin' (id da bacia por frame, 0 = fora).
        """
        data = self.combined_data()
        result = self.calculate_free_energy(data)
        G, X, Y = result['G_original'], result['X_original'], result['Y_original']

        if self.max_energy is not None and threshold > self.max_energy:
            raise ValueError(
                f"--energy ({threshold}) is above --max_energy ({self.max_energy}); the surface "
                f"is capped there, so every capped cell would merge into one basin. "
                f"Raise --max_energy or lower --energy."
            )

        structure = ndimage.generate_binary_structure(2, connectivity)
        below = np.less_equal(G, threshold, where=np.isfinite(G),
                              out=np.zeros(G.shape, dtype=bool))
        raw_labels, n_found = ndimage.label(below, structure=structure)

        energies = self.frame_energies()
        i, j, inside = self._grid_indices(result)
        raw_frame_label = np.where(inside, raw_labels[i, j], 0)
        saddles = self._basin_saddles(G, raw_labels)

        candidates = []
        for lab in range(1, n_found + 1):
            cells = raw_labels == lab
            members = np.flatnonzero(raw_frame_label == lab)
            if members.size < min_frames:
                continue
            gi, gj = np.unravel_index(np.nanargmin(np.where(cells, G, np.nan)), G.shape)
            rep = members[np.argmin(energies[members])]
            neighbours = [h for (a, b), h in saddles.items() if lab in (a, b)]
            candidates.append({
                'raw_label': lab,
                'G_min': float(G[gi, gj]),
                'depth': float(min(neighbours) - G[gi, gj]) if neighbours else float('inf'),
                'cv1_min': float(X[gi, gj]),
                'cv2_min': float(Y[gi, gj]),
                'n_frames': int(members.size),
                'population': float(members.size / len(energies)),
                'rep_frame': int(self.frames[rep]),
                'rep_cv1': float(self.proj1_data_original[rep]),
                'rep_cv2': float(self.proj2_data_original[rep]),
                'rep_energy': float(energies[rep]),
                'area_cells': int(cells.sum()),
            })

        # Bacia 1 e sempre a de menor energia
        candidates.sort(key=lambda b: b['G_min'])

        labels = np.zeros_like(raw_labels)
        frame_basin = np.zeros_like(raw_frame_label)
        for new_id, basin in enumerate(candidates, start=1):
            labels[raw_labels == basin['raw_label']] = new_id
            frame_basin[raw_frame_label == basin['raw_label']] = new_id
            basin['id'] = new_id
            del basin['raw_label']

        return {'basins': candidates, 'labels': labels, 'frame_basin': frame_basin,
                'threshold': threshold, 'method': 'connected'}

    def print_basins(self, basins):
        """Resumo das bacias no terminal."""
        if not basins:
            print("No basin found below the energy threshold.")
            return
        print(f"\n{len(basins)} basin(s) identified:\n")
        head = (f"{'#':>3} {'CV1':>10} {'CV2':>10} {'dG':>8} {'depth':>8} "
                f"{'frames':>8} {'pop%':>7} {'rep frame':>10}")
        print(head)
        print("-" * len(head))
        for b in basins:
            depth = b.get('depth', float('inf'))
            depth_txt = "  -" if not np.isfinite(depth) else f"{depth:>8.2f}"
            print(f"{b['id']:>3} {b['cv1_min']:>10.3f} {b['cv2_min']:>10.3f} "
                  f"{b['G_min']:>8.2f} {depth_txt:>8} {b['n_frames']:>8d} "
                  f"{100 * b['population']:>7.2f} {b['rep_frame']:>10d}")
        print("\ndG    = free energy of the basin minimum, relative to the global minimum")
        print("depth = height of the lowest saddle to a neighbouring basin, above that minimum\n")

    @staticmethod
    def print_persistence(result):
        """Mostra o que --basin_min_depth fundiu, para que a escolha nao fique invisivel."""
        spectrum = result.get('persistence')
        if not spectrum:
            return
        cut = result['min_depth']
        merged = [d for d in spectrum if d < cut]
        how = "auto" if result.get('auto_depth') else "set"
        print(f"Watershed found {result['n_raw_minima']} local minima; "
              f"--basin_min_depth {cut:.2f} kJ/mol ({how}) merged {len(merged)} of them.")
        shown = ", ".join(f"{d:.2f}" for d in spectrum[:12])
        print(f"Barrier depths (deepest first): {shown}"
              + (" ..." if len(spectrum) > 12 else ""))
        # Um degrau largo no espectro separa estrutura real de ruido do KDE
        if len(spectrum) > 1:
            gaps = [(spectrum[i] - spectrum[i + 1], i) for i in range(len(spectrum) - 1)]
            gap, i = max(gaps)
            if gap > 0.2:
                print(f"Largest step is {spectrum[i]:.2f} -> {spectrum[i + 1]:.2f} kJ/mol: "
                      f"--basin_min_depth between them keeps {i + 2} basins.")
        print()

    def save_basins(self, basins, filename='basins.tsv'):
        """Grava a tabela de bacias com o frame representativo de cada uma."""
        columns = ['cluster', 'cv1_min', 'cv2_min', 'G_min_kJ_mol', 'depth_kJ_mol',
                   'n_frames', 'population_fraction', 'rep_frame', 'rep_cv1', 'rep_cv2',
                   'rep_energy_kJ_mol', 'area_cells']
        rows = [[b['id'], b['cv1_min'], b['cv2_min'], b['G_min'],
                 b.get('depth', float('inf')), b['n_frames'],
                 b['population'], b['rep_frame'], b['rep_cv1'], b['rep_cv2'],
                 b['rep_energy'], b['area_cells']] for b in basins]
        np.savetxt(filename,
                   np.array(rows, dtype=np.float64) if rows else np.empty((0, len(columns))),
                   delimiter='\t',
                   fmt=['%d', '%.6f', '%.6f', '%.6f', '%.6f', '%d', '%.6f',
                        '%d', '%.6f', '%.6f', '%.6f', '%d'],
                   header='\t'.join(columns),
                   comments='')
        print(f"Basin summary saved in '{filename}'.")

    def _view_limits(self, result, pad=0.03):
        """Limites de eixo ajustados a regiao amostrada, respeitando --xlim/--ylim."""
        X, Y, mask = result['X_original'], result['Y_original'], result['mask']
        xs, ys = X[mask], Y[mask]
        span_x = (xs.max() - xs.min()) * pad
        span_y = (ys.max() - ys.min()) * pad
        x_lo = self.xlim_inf if self.xlim_inf is not None else xs.min() - span_x
        x_hi = self.xlim_sup if self.xlim_sup is not None else xs.max() + span_x
        y_lo = self.ylim_inf if self.ylim_inf is not None else ys.min() - span_y
        y_hi = self.ylim_sup if self.ylim_sup is not None else ys.max() + span_y
        return (x_lo, x_hi), (y_lo, y_hi)

    def _surface_3d(self, ax, result):
        """Superficie 3D restrita a regiao amostrada.

        Cada quadrilatero so e desenhado quando os seus quatro vertices estao dentro da
        mascara; os demais ficam transparentes. Isso evita tanto os buracos que o NaN
        produz na triangulacao quanto uma tampa opaca sobre a paisagem.
        """
        X, Y = result['X_original'], result['Y_original']
        G, mask = result['G_original'], result['mask']

        ceiling = float(np.nanmax(G))
        # Achatar no teto e so entao recortar. Cortando exatamente na borda da mascara,
        # o corte cai no meio da parede: como a regiao amostrada e uma faixa diagonal e a
        # grade e ortogonal, cada linha termina numa altura diferente e a silhueta vira
        # uma serra de barbatanas. Estendendo o recorte alguns passos alem da mascara, o
        # corte acontece em terreno ja plano e a borda fica limpa.
        Z = np.minimum(result['G_unmasked'], ceiling)
        Z = np.minimum(ndimage.gaussian_filter(Z, sigma=1.0), ceiling)
        render = ndimage.binary_dilation(mask, iterations=6)

        norm = plt.Normalize(0.0, ceiling)
        cmap = self.get_cmap()
        colors = cmap(norm(Z))
        quad = render[:-1, :-1] & render[1:, :-1] & render[:-1, 1:] & render[1:, 1:]
        colors[:-1, :-1, 3] = quad.astype(float)

        surf = ax.plot_surface(X, Y, Z, facecolors=colors, shade=False,
                               rstride=1, cstride=1, linewidth=0, antialiased=True)
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        mappable.set_array([])
        return surf, mappable, ceiling

    def _annotate_basins(self, ax, basins, three_d=False):
        """Marca o minimo de cada bacia com o seu numero e devolve os itens da legenda.

        O rotulo sai deslocado do marcador, com linha de chamada, para que bacias vizinhas
        nao sobreponham os numeros.
        """
        handles = []
        # Deslocamentos alternados evitam colisao entre rotulos de bacias proximas
        offsets = [(16, 14), (-20, 14), (16, -18), (-20, -18), (22, 0), (-26, 0)]

        for b in basins:
            color = self.discreet_colors[(b['id'] - 1) % len(self.discreet_colors)]
            marker = self.discreet_markers[(b['id'] - 1) % len(self.discreet_markers)]

            if three_d:
                # Um marcador desenhado em 3D some atras da propria superficie quando a
                # bacia fica do lado oposto ao da camera, e zorder nao ajuda porque o
                # mplot3d ordena por profundidade. Projetando o minimo para o plano da
                # figura e desenhando ali, tanto a forma quanto o rotulo ficam sempre
                # por cima, com uma linha de chamada apontando o ponto — igual ao 2D.
                x2, y2, _ = proj3d.proj_transform(b['cv1_min'], b['cv2_min'],
                                                  b['G_min'], ax.get_proj())
                anchor = (x2, y2)
                ax.add_artist(Line2D([x2], [y2], marker=marker, color=color,
                                     markersize=11, markeredgecolor='white',
                                     markeredgewidth=1.4, linestyle='none',
                                     transform=ax.transData, zorder=9))
            else:
                ax.scatter([b['cv1_min']], [b['cv2_min']], color=color, marker=marker,
                           s=130, edgecolors='white', linewidths=1.4, zorder=6)
                anchor = (b['cv1_min'], b['cv2_min'])

            ax.annotate(
                str(b['id']),
                xy=anchor,
                xytext=offsets[(b['id'] - 1) % len(offsets)],
                textcoords='offset points',
                fontsize=10, fontweight='bold', color='black',
                ha='center', va='center', zorder=10,
                bbox=dict(boxstyle='circle,pad=0.25', facecolor='white',
                          edgecolor=color, linewidth=1.4),
                arrowprops=dict(arrowstyle='-', color=color, linewidth=1.2,
                                shrinkA=0, shrinkB=6),
            )

            depth = b.get('depth', float('inf'))
            depth_txt = "" if not np.isfinite(depth) else f", depth {depth:.1f}"
            handles.append(plt.Line2D([0], [0], marker=marker, color='none',
                                      markerfacecolor=color, markeredgecolor='white',
                                      markeredgewidth=1.2, markersize=9,
                                      label=(f"{b['id']}  dG {b['G_min']:.2f} kJ/mol"
                                             f"{depth_txt}  |  "
                                             f"{100 * b['population']:.1f}% of frames  |  "
                                             f"frame {b['rep_frame']}")))
        return handles

    def _profile_1d(self, data):
        """Perfil de energia livre de uma CV, a partir do histograma 1D.

        Bins vazios ficam NaN em vez de receberem um piso artificial de densidade.
        Piso de 1e-10 virava um pico de ~57 kJ/mol no grafico, indistinguivel de uma
        barreira medida quando na verdade significa "nenhum frame aqui".
        """
        hist, bin_edges = np.histogram(data, bins=self.bins, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        with np.errstate(divide='ignore'):
            G = -self.kB * self.temperature * np.log(np.where(hist > 0, hist, np.nan))
        return bin_centers, G - np.nanmin(G)

    def boltzmann_inversion(self, data_list, titles, threshold=None):

        fig_combined, axs_combined = plt.subplots(1, len(data_list), figsize=(14, 5),
                                                  sharey=True)
        # Inicializa a lista de elementos da legenda com o elemento para a linha de energia livre
        legend_elements = [plt.Line2D([0], [0], color='red', lw=2, label='Free energy')]

        for idx, (ax, data, title) in enumerate(zip(axs_combined, data_list, titles)):
            bin_centers, G_min_normalized = self._profile_1d(data)

            # Plotar a linha de energia livre
            ax.plot(bin_centers, G_min_normalized,
                    color='red', label='Free energy'
                    if idx == 0 else "_nolegend_"
                    )

            ax.set_xlabel(title)
            ax.grid(True, alpha=0.25, linewidth=0.5)

        if threshold is not None and hasattr(self, 'discrete') and self.discrete:
            discrete_intervals = np.arange(0, threshold, self.discrete)
            for i, interval in enumerate(discrete_intervals):
                end = interval + self.discrete if i < len(discrete_intervals) - 1 else threshold
                # Evita adicionar intervalos redundantes na legenda
                if end > interval:
                    legend_elements.append(plt.Line2D([0], [0], 
                                                      marker=self.discreet_markers[i % len(self.discreet_markers)], 
                                                      color='w', markerfacecolor=self.discreet_colors[i % len(self.discreet_colors)], 
                                                      markersize=10, label=f'{interval:.1f}-{end:.1f} kJ/mol'))

                # Plotar os pontos para cada intervalo nos gráficos
                for ax, data in zip(axs_combined, data_list):
                    bin_centers, G_min_normalized = self._profile_1d(data)
                    mask = (G_min_normalized >= interval) & (G_min_normalized < end)
                    if np.any(mask):
                        ax.scatter(bin_centers[mask], G_min_normalized[mask], 
                                   color=self.discreet_colors[i % len(self.discreet_colors)], 
                                   marker=self.discreet_markers[i % len(self.discreet_markers)], 
                                   s=50)

        # Titulo de legenda so faz sentido quando ha intervalos a listar
        axs_combined[-1].legend(handles=legend_elements,
                                loc='upper right', frameon=False,
                                title="Energy intervals" if len(legend_elements) > 1 else None)

        axs_combined[0].set_ylabel('Free Energy (kJ/mol)')
        plt.suptitle('Free Energy Profile of Each Collective Variable')
        plt.tight_layout()
        plt.savefig('Combined_Free_Energy_Profile_Normalized.png',
                    dpi=self.dpi, bbox_inches='tight')
        plt.show()
        

    def plot_histogram(self, data_list, titles):
        plt.figure(figsize=(7 * len(data_list), 5))
        
        # Contagens em porcentagem do bin mais alto entre as duas CVs. O eixo x fica nas
        # unidades originais da CV, para que este grafico possa ser lido lado a lado com
        # o perfil de energia livre, que usa os mesmos valores.
        all_counts = [np.histogram(data, bins=self.bins)[0] for data in data_list]
        total_counts_max = max(counts.max() for counts in all_counts)

        for i, (data, title) in enumerate(zip(data_list, titles)):
            counts, bin_edges = np.histogram(data, bins=self.bins)
            normalized_counts = (counts / total_counts_max) * 100
            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

            ax = plt.subplot(1, len(data_list), i + 1)
            ax.bar(bin_centers,
                   normalized_counts,
                   width=(bin_edges[1] - bin_edges[0]),
                   alpha=0.85,
                   color='#3a7d44'
                   )

            ax.set_ylim(0, 100)  # Definir explicitamente o eixo Y de 0 a 100%
            ax.set_yticks(np.arange(0, 101, 20))
            ax.set_xlabel(title)
            if i == 0:
                ax.set_ylabel('Frequency (% of tallest bin)')
            ax.grid(True, alpha=0.25, linewidth=0.5)

        plt.suptitle('Distribution of Each Collective Variable')
        plt.tight_layout()
        plt.savefig('histograms_normalized_side_by_side.png',
                    dpi=self.dpi, bbox_inches='tight')
        plt.show()


    def cv_by_frame(self, data_list, titles):
        frames = np.arange(len(data_list[0]))  # Assumindo que todos os conjuntos têm o mesmo número de frames
        for data, title in zip(data_list, titles):
            plt.figure(figsize=(10, 6))
            plt.plot(frames, data, label=title)
            plt.xlabel('Frame')
            plt.ylabel(title)
            plt.title(f'CV by Frame - {title}')
            plt.legend()
            # plt.savefig(f'cv_by_frame_absolute_{title.replace(" ", "_")}.png')
            plt.close()

        # Plot combinado dos valores relativos, ajustado para exibir de 0 a 100
        plt.figure(figsize=(10, 6))
        for data, title in zip(data_list, titles):
            # Normalizando os dados de 0 a 1 e multiplicando por 100 para ajustar a escala de 0 a 100
            data_normalized = ((data - np.min(data)) / (np.max(data) - np.min(data))) * 100
            plt.plot(frames, data_normalized, label=title)
        plt.xlabel('Frame')
        plt.ylabel('CV (%)')  # Atualização da etiqueta do eixo Y para refletir a escala de porcentagem
        plt.title('CV by Frame - Combined Normalized')
        plt.legend()
        plt.ylim(0, 100)  # Garantindo que o eixo Y esteja limitado de 0 a 100
        plt.savefig('cv_by_frame_combined_normalized.png',
                    dpi=self.dpi, bbox_inches='tight')
        plt.show()


    def plot_energy_landscape(self, threshold, titles=['CV1', 'CV2'], basins=None,
                              xlim_inf=None, xlim_sup=None, ylim_inf=None, ylim_sup=None):

        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)
        X, Y, G = result['X_original'], result['Y_original'], result['G_original']

        cmap = self.get_cmap()
        fig, ax = plt.subplots(figsize=(9, 6.5))
        # O fundo recebe a cor do topo da paleta, continuando a escala em vez de cortar
        # a figura com um vazio branco. Fora da regiao amostrada a energia livre e de
        # fato mais alta (foi por isso que nada foi amostrado ali); o que nao se sabe e
        # quanto. A fronteira da amostragem vai marcada por uma linha tracejada.
        ax.set_facecolor(cmap(1.0))

        # Niveis so na faixa util. Passos redondos, cerca de uma duzia de bandas.
        G_max = float(np.nanmax(G))
        step = max(round(G_max / 12, 1), 0.1)
        levels = np.arange(0, G_max + step, step)

        cont = ax.contourf(X, Y, G, levels=levels, cmap=cmap, extend='max')
        # Linhas finas e claras: dao relevo sem virar emaranhado preto
        ax.contour(X, Y, G, levels=levels, colors='white', linewidths=0.4, alpha=0.45)
        # Ate onde a amostragem chega. Discreto de proposito: uma linha forte aqui
        # recorta a faixa amostrada do fundo e a figura passa a parecer colada.
        ax.contour(X, Y, result['mask'].astype(float), levels=[0.5],
                   colors='white', linewidths=0.7, linestyles='dashed', alpha=0.3)

        handles = []
        legend_title = None

        if basins:
            # Cristas do watershed, apenas entre duas bacias reais. Fora delas o campo
            # de rotulos vai a NaN para que nenhuma linha seja tracada na borda externa:
            # contorno da bacia, borda da mascara e limite dos rotulos cairiam quase uns
            # sobre os outros e formariam um halo branco grosso, que e o que faz a faixa
            # amostrada parecer um adesivo colado sobre o fundo.
            n_basins = len(basins['basins'])
            if n_basins > 1:
                ridges = np.where(basins['labels'] > 0,
                                  basins['labels'].astype(float), np.nan)
                ax.contour(X, Y, ridges, levels=np.arange(1.5, n_basins + 0.5),
                           colors='white', linewidths=1.0, alpha=0.65, zorder=4)
            handles = self._annotate_basins(ax, basins['basins'])
            legend_title = (f"Basins ({basins.get('method', 'watershed')})"
                            + (f", minimum below {basins['threshold']:g} kJ/mol"
                               if basins.get('threshold') is not None else ""))

        elif self.discrete is not None and threshold is not None:
            # Modo legado: marca nos da grade por faixa de energia (nao sao conformacoes)
            discrete_intervals = np.arange(0, threshold, self.discrete)

            for i, interval in enumerate(discrete_intervals):
                end = min(interval + self.discrete, threshold)
                mask = (result['G_original'].flatten() <= end) & (result['G_original'].flatten() > interval)
                X_flat, Y_flat = result['X_original'].flatten(), result['Y_original'].flatten()

                if np.any(mask):
                    ax.scatter(X_flat[mask],
                               Y_flat[mask],
                               color=self.discreet_colors[i % len(self.discreet_colors)],
                               marker=self.discreet_markers[i % len(self.discreet_markers)],
                               label=f'{interval:.1f}-{end:.1f} KJ/mol'
                               )
            handles, _ = ax.get_legend_handles_labels()
            legend_title = "Energy intervals"

        # Colorbar colada ao eixo; legenda abaixo, para nenhuma das duas roubar area do plot
        cbar = fig.colorbar(cont, ax=ax, pad=0.02, fraction=0.046)
        cbar.set_label('Free energy (kJ/mol)')
        cbar.outline.set_linewidth(0.5)

        (x_lo, x_hi), (y_lo, y_hi) = self._view_limits(result)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_xlabel(titles[0])
        ax.set_ylabel(titles[1])
        ax.set_title('Free Energy Landscape', pad=10)

        if handles:
            ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.13),
                      borderaxespad=0., title=legend_title, fontsize=8,
                      title_fontsize=9, frameon=False,
                      ncol=1 if len(handles) <= 4 else 2)

        fig.savefig('Free_energy_landscape_with_discrete_points.png',
                    dpi=self.dpi, bbox_inches='tight')
        plt.show()

    def plot_3D_energy_landscape(self, threshold=None, titles=['CV1', 'CV2'], basins=None,
                                 xlim_inf=None, xlim_sup=None, ylim_inf=None, ylim_sup=None):
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)

        fig = plt.figure(figsize=(10, 7.5))
        ax = fig.add_subplot(111, projection='3d')

        surf, mappable, ceiling = self._surface_3d(ax, result)
        ax.set_zlim(0, ceiling)

        ax.set_xlabel(titles[0], labelpad=10)
        ax.set_ylabel(titles[1], labelpad=10)
        ax.set_zlabel('Free energy (kJ/mol)', labelpad=8)
        ax.set_title('3D Free Energy Landscape', pad=14)
        ax.view_init(elev=32, azim=-125)
        ax.set_box_aspect((1, 1, 0.55))
        # Painel de fundo discreto: a superficie e o objeto de interesse
        for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
            pane.pane.set_alpha(0.04)
            pane.pane.set_edgecolor('none')

        # Limites antes das anotacoes: a projecao 3D -> 2D usada para posicionar os
        # rotulos depende deles e da camera.
        (x_lo, x_hi), (y_lo, y_hi) = self._view_limits(result)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)

        handles = []
        legend_title = None

        if basins:
            handles = self._annotate_basins(ax, basins['basins'], three_d=True)
            legend_title = (f"Basins ({basins.get('method', 'watershed')})"
                            + (f", minimum below {basins['threshold']:g} kJ/mol"
                               if basins.get('threshold') is not None else ""))

        elif self.discrete is not None and threshold is not None:
            # Modo legado: marca nos da grade por faixa de energia (nao sao conformacoes)
            discrete_intervals = np.arange(0, threshold, self.discrete)
            for i, interval in enumerate(discrete_intervals):
                end = min(interval + self.discrete, threshold)
                mask = (result['G_original'].flatten() <= end) & (result['G_original'].flatten() > interval)
                X_flat, Y_flat, Z_flat = result['X_original'].flatten(), result['Y_original'].flatten(), result['G_original'].flatten()

                if np.any(mask):
                    ax.scatter(X_flat[mask],
                               Y_flat[mask],
                               Z_flat[mask],
                               color=self.discreet_colors[i % len(self.discreet_colors)],
                               marker=self.discreet_markers[i % len(self.discreet_markers)],
                               label=f'{interval:.1f}-{end:.1f} KJ/mol')
            handles, _ = ax.get_legend_handles_labels()
            legend_title = "Energy intervals"

        cbar = fig.colorbar(mappable, ax=ax,
                            shrink=0.55,
                            aspect=18,
                            pad=0.10,
                            label='Free energy (kJ/mol)'
                            )
        cbar.outline.set_linewidth(0.5)

        # Evita o UserWarning quando nao ha artistas rotulados
        if handles:
            ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.02),
                      title=legend_title, fontsize=8, title_fontsize=9, frameon=False,
                      ncol=1 if len(handles) <= 4 else 2)

        fig.savefig('3D_landscape.png', dpi=self.dpi, bbox_inches='tight')
        plt.show()


    def plot_threshold_points(self, ax, result, lower_bound, upper_bound, color, label):
        G_flat = result['G_original'].flatten()
        energy_mask = (G_flat >= lower_bound) & (G_flat < upper_bound)

        if any(energy_mask):
            X_flat, Y_flat = result['X_original'].flatten(), result['Y_original'].flatten()
            ax.scatter(X_flat[energy_mask], Y_flat[energy_mask], G_flat[energy_mask], color=color, s=20, label=label)


    def create_3D_gif(self, gif_filename='energy_landscape_3D.gif', n_angles=10, elevation=15, duration_per_frame=0.01, titles=['CV1', 'CV2'], xlim_inf=None, xlim_sup=None, ylim_inf=None, ylim_sup=None):
        
        temp_dir = tempfile.mkdtemp()  # Cria um diretório temporário para armazenar os frames
        filenames = []

        # Utiliza a função calculate_free_energy para obter os dados
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)

        if n_angles < 1:
            raise ValueError(f"--gif_angles must be at least 1, got {n_angles}.")
        # Angulos igualmente espacados; nao usar passo inteiro, que quebra para n > 360
        angles = np.linspace(0, 360, int(n_angles), endpoint=False)

        (x_lo, x_hi), (y_lo, y_hi) = self._view_limits(result)

        for i, angle in enumerate(angles):
            fig = plt.figure(figsize=(10, 7.5))
            ax = fig.add_subplot(111, projection='3d')

            surf, mappable, vmax = self._surface_3d(ax, result)
            vmin = 0.0

            ax.view_init(elev=elevation, azim=angle)
            ax.set_xlabel(titles[0], labelpad=10)
            ax.set_ylabel(titles[1], labelpad=10)
            ax.set_zlabel('Free energy (kJ/mol)', labelpad=8)
            ax.set_title('3D Free Energy Landscape', pad=14)
            ax.set_xlim(x_lo, x_hi)
            ax.set_ylim(y_lo, y_hi)
            ax.set_zlim(vmin, vmax)
            ax.set_box_aspect((1, 1, 0.55))
            for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
                pane.pane.set_alpha(0.04)
                pane.pane.set_edgecolor('none')

            # Colorbar em todos os frames: aparecer so em alguns faz o GIF "piscar"
            fig.colorbar(mappable, ax=ax, shrink=0.55, aspect=18, pad=0.10,
                         label='Free energy (kJ/mol)')

            frame_filename = os.path.join(temp_dir, f"frame_{i:03d}.png")
            fig.savefig(frame_filename, dpi=110)
            filenames.append(frame_filename)
            plt.close(fig)

        # imageio >= 2.28 interpreta duration de GIF em milissegundos
        with imageio.get_writer(gif_filename, mode='I',
                                duration=max(duration_per_frame * 1000, 20),
                                loop=0) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)

        shutil.rmtree(temp_dir)  # Limpa os arquivos temporários

        # Abrir o GIF gerado automaticamente
        self.open_gif(gif_filename)


    def open_gif(self, gif_filename):
        if platform.system() == 'Windows':
            os.startfile(gif_filename)
        elif platform.system() == 'Darwin':  # macOS
            subprocess.run(['open', gif_filename])
        else:  # Assume Linux ou outra plataforma Unix-like
            subprocess.run(['xdg-open', gif_filename])

    @staticmethod
    def help():
        help_text = """
        Usage:
            free_energy_landscape path/to/cv1_data.txt path/to/cv2_data.txt

        Optional arguments:
            --temperature           [int]       Simulation temperature in Kelvin (default: 300K)
            --kb                    [float]     Boltzmann constant in kJ/(mol·K) (default: 8.314e-3)
            --energy                [float]     Energy threshold (KJ/mol). With the watershed
                                                method it only filters which basins are reported
                                                (those whose minimum lies below it) (default: None)
            --basin_method          [str]       'watershed' (default) assigns every frame to the
                                                minimum it descends to; 'connected' groups only
                                                the cells below --energy
            --basin_min_depth       [float]     Merge basins separated by a barrier shallower than
                                                this (KJ/mol). Default: kB*T, the thermal scale
            --basin_connectivity    [int]       'connected' only. 1 = 4 neighbours, 2 = 8 (default: 2)
            --basin_min_frames      [int]       Discard basins holding fewer sampled frames than
                                                this, as KDE artefacts (default: 1)
            --mask_min_count        [float]     A grid cell counts as sampled when the KDE expects
                                                at least this many frames inside it; the rest is
                                                left blank instead of extrapolated (default: 1)
            --cmap                  [str]       Matplotlib colormap, or 'classic' for the 1.x
                                                blue-to-red map (default: viridis_r)
            --dpi                   [int]       Resolution of the saved figures (default: 200)
            --discretize            [float]     LEGACY. Colours individual grid nodes by energy
                                                band instead of grouping them into basins. Those
                                                nodes are grid cells, not conformations (default: None)
            --bins_energy_histogram [int]       Bins for energy histogram (default: 100)
            --kde_bandwidth         [float]     Bandwidth for kernel density estimation (default: None)
            --max_energy            [float]     Upper cap (KJ/mol) for the plotted energy surface
                                                (default: none). Only affects the colour scale of
                                                the plots, never the .tsv values.
            --grid                  [int]       Grid resolution per axis for the landscape (default: 100)
            --n_jobs                [int]       Parallel workers for KDE evaluation (default: all CPUs).
                                                Results are identical for any value.
            --names                 [str] [str] Names for the collective variables (default: CV1, CV2)
            --gif_angles            [int]       Angles for 3D GIF rotation (default: 10)
            --gif_elevation         [int]       Elevation angle for the 3D GIF (default: 10)
            --gif_duration          [float]     Duration per frame in the GIF in seconds (default: 0.1)
            --xlim_inf              [float]     Lower limit for the x-axis (default: None)
            --xlim_sup              [float]     Upper limit for the x-axis (default: None)
            --ylim_inf              [float]     Lower limit for the y-axis (default: None)
            --ylim_sup              [float]     Upper limit for the y-axis (default: None)

        Example:
            freeEnergyLandscape.py proj1Out.txt proj2Out.txt --energy 3.0 --names "CV1 (Angle)" "CV2 (Distance)"

        With --energy, the region below the threshold is grouped into basins (connected
        components of the free energy surface). Each basin is numbered on the plots and
        reported in 'basins.tsv' together with its representative frame: the sampled
        conformation with the lowest free energy inside that basin.

        """
        print(help_text)

    def calculate_and_save_free_energy(self, threshold=None, filename='discrete_values_energy_frames.tsv'):
        """Energia livre por frame amostrado, usando o MESMO kernel global das superficies."""

        # Verifica se os dados foram carregados
        if self.proj1_data_original is None or self.proj2_data_original is None:
            raise ValueError("Data not loaded. Run load_data first.")

        frames = self.frames

        # Prepara os dados combinados
        combined_data = np.vstack((self.proj1_data_original, self.proj2_data_original)).T

        # Um unico kernel ajustado sobre toda a trajetoria; a paralelizacao divide
        # apenas os pontos de avaliacao, portanto o resultado independe do numero de CPUs.
        # (frame_energies reutiliza esse mesmo kernel.)
        G_normalized = self.frame_energies()

        # Se as bacias ja foram identificadas, cada frame carrega o id da sua bacia
        basin_of_frame = self._basins['frame_basin'] if self._basins else None

        # Aplica o threshold, se especificado
        if threshold is not None:
            keep = G_normalized <= threshold
        else:
            keep = np.ones(G_normalized.shape, dtype=bool)

        columns = [frames[keep], self.proj1_data_original[keep],
                   self.proj2_data_original[keep], G_normalized[keep]]
        header = ['frame', 'cv1', 'cv2', 'energy']
        fmt = ['%d', '%.6f', '%.6f', '%.6f']
        if basin_of_frame is not None:
            columns.append(basin_of_frame[keep])
            header.append('cluster')
            fmt.append('%d')

        # Prepara os dados para salvamento
        data_to_save = np.column_stack(columns)

        # Ordena os dados pela energia
        data_to_save = data_to_save[data_to_save[:, 3].argsort()]

        # Salva os dados em um arquivo .tsv
        np.savetxt(filename,
                   data_to_save,
                   delimiter='\t',
                   fmt=fmt,
                   header='\t'.join(header),
                   comments=''
                   )

        print(f"Energy data saved in '{filename}'.")


    def main(self, energy_threshold, cv_names, n_angles, elevation, duration_per_frame):

        print("Loading data...\n")
        self.load_data()

        print("Data loaded successfully!")
        print(f"CV1: {self.cv1_path}, CV2: {self.cv2_path}\n")

        print("Plotting free energy profiles...\n")
        self.boltzmann_inversion(
            data_list=[self.proj1_data_original, self.proj2_data_original],
            titles=cv_names,
            threshold=energy_threshold
            )

        print("Plotting histograms...\n")
        self.plot_histogram(
            data_list=[self.proj1_data_original, self.proj2_data_original],
            titles=cv_names
            )
        print("Plotting Collective variables normalized by frame...\n")
        self.cv_by_frame(
            data_list=[self.proj1_data_original, self.proj2_data_original],
            titles=cv_names
            )
        
        print("Successfully generated, free energy profiles, histograms and normalized collective variables.\n")

        basins = None
        if self.discrete is None and (energy_threshold is not None
                                      or self.basin_method == 'watershed'):
            print(f"Identifying basins ({self.basin_method})...\n")
            basins = self.identify_basins(energy_threshold,
                                          connectivity=self.basin_connectivity,
                                          min_frames=self.basin_min_frames)
            self._basins = basins
            self.print_persistence(basins)
            self.print_basins(basins['basins'])
            self.save_basins(basins['basins'])

        print("Plotting the free energy landscape...\n")
        self.plot_energy_landscape(
            threshold=energy_threshold, titles=cv_names, basins=basins
            )
        print("Plot successfully generated.\n")

        print("Plotting the free energy landscape in 3D...\n")
        self.plot_3D_energy_landscape(
            threshold=energy_threshold, titles=cv_names, basins=basins
                                      )
        print("Plotting 3D gif...\n")

        self.create_3D_gif(
            n_angles=n_angles, elevation=elevation,
            duration_per_frame=duration_per_frame,
            titles=cv_names
                           )
        print("3D plot successfully generated.\n")

        # Após o uso final dos dados, limpe-os para liberar memória
        self.cached_results = None
        self._cache_key = None

def main():
    # Definindo valores padrão
    t = 300                     # --temperature           [int] [Kelvin]
    kB = 8.314e-3               # --kb                    [float] [kJ/(mol·K)]
    energy_threshold = None     # --energy                [float] [kJ/mol]
    bins_energy_histogram = 100 # --bins_energy_histogram [int]
    kde_bandwidth_cv = None     # --kde_bandwidth         [float]
    cv_names = ['CV1', 'CV2']   # --name                  [str] [str]
    n_angles = 10               # --gif_angles            [int]
    elevation = 10              # --gif_elevation         [int]
    duration_per_frame = 0.1    # --gif_duration          [float]
    discrete_val = None         # --discrete              [float]
    max_energy = None           # --max_energy            [float] [kJ/mol] ou None
    grid_size = 100             # --grid                  [int]
    n_jobs = None               # --n_jobs                [int]
    basin_method = 'watershed'  # --basin_method          [str]
    basin_connectivity = 2      # --basin_connectivity    [int]
    basin_min_frames = 1        # --basin_min_frames      [int]
    basin_min_depth = None      # --basin_min_depth       [float] [kJ/mol]
    mask_min_count = 1.0        # --mask_min_count        [float]
    cmap = 'viridis_r'          # --cmap                  [str]
    dpi = 200                   # --dpi                   [int]
    xlim_inf = xlim_sup = ylim_inf = ylim_sup = None  # Inicialização padrão

    if len(sys.argv) >= 2 and sys.argv[1] in ("-h", "--help"):
        FreeEnergyLandscape.help()
        sys.exit(0)

    if len(sys.argv) >= 3:
        cv1_path, cv2_path = sys.argv[1], sys.argv[2]

        # Processar argumentos adicionais como pares chave-valor
        i = 3
        while i < len(sys.argv):
            key = sys.argv[i]
            needed = 3 if key == "--names" else 2
            if key.startswith("--") and i + needed > len(sys.argv):
                print(f"Option {key} requires {needed - 1} value(s).")
                sys.exit(1)
            if key == "--temperature":
                t = float(sys.argv[i + 1])
                i += 2
            elif key == "--kb":
                kB = float(sys.argv[i + 1])
                i += 2
            elif key == "--energy":
                energy_threshold = float(sys.argv[i + 1])
                i += 2
            elif key == "--discretize":
                discrete_val = float(sys.argv[i + 1])  
                i += 2
            elif key == "--bins_energy_histogram":
                bins_energy_histogram = int(sys.argv[i + 1])
                i += 2
            elif key == "--kde_bandwidth":
                kde_bandwidth_cv = float(sys.argv[i + 1]) if sys.argv[i + 1].lower() != "none" else None
                i += 2
            elif key == "--max_energy":
                max_energy = float(sys.argv[i + 1]) if sys.argv[i + 1].lower() != "none" else None
                i += 2
            elif key == "--grid":
                grid_size = int(sys.argv[i + 1])
                i += 2
            elif key == "--n_jobs":
                n_jobs = int(sys.argv[i + 1])
                i += 2
            elif key == "--basin_method":
                basin_method = sys.argv[i + 1]
                i += 2
            elif key == "--basin_connectivity":
                basin_connectivity = int(sys.argv[i + 1])
                i += 2
            elif key == "--basin_min_frames":
                basin_min_frames = int(sys.argv[i + 1])
                i += 2
            elif key == "--basin_min_depth":
                basin_min_depth = float(sys.argv[i + 1])
                i += 2
            elif key == "--mask_min_count":
                mask_min_count = float(sys.argv[i + 1])
                i += 2
            elif key == "--cmap":
                cmap = sys.argv[i + 1]
                i += 2
            elif key == "--dpi":
                dpi = int(sys.argv[i + 1])
                i += 2
            elif key == "--names":
                cv_names = [sys.argv[i + 1], sys.argv[i + 2]]
                i += 3
            elif key == "--gif_angles":
                n_angles = int(sys.argv[i + 1])
                i += 2
            elif key == "--gif_elevation":
                elevation = int(sys.argv[i + 1])
                i += 2
            elif key == "--gif_duration":
                duration_per_frame = float(sys.argv[i + 1])
                i += 2

            elif key == "--xlim_inf":
                xlim_inf = float(sys.argv[i + 1])
                i += 2
            elif key == "--xlim_sup":
                xlim_sup = float(sys.argv[i + 1])
                i += 2
            elif key == "--ylim_inf":
                ylim_inf = float(sys.argv[i + 1])
                i += 2
            elif key == "--ylim_sup":
                ylim_sup = float(sys.argv[i + 1])
                i += 2

            else:
                print(f"Unrecognized option: {key}")
                sys.exit(1)
    else:
        FreeEnergyLandscape.help()
        sys.exit(1)

    try:
        fel = FreeEnergyLandscape(cv1_path, cv2_path, t, kB, 
                                bins=bins_energy_histogram, 
                                kde_bandwidth=kde_bandwidth_cv, 
                                cv_names=cv_names, 
                                discrete=discrete_val,
                                xlim_inf=xlim_inf, xlim_sup=xlim_sup,
                                ylim_inf=ylim_inf, ylim_sup=ylim_sup,
                                max_energy=max_energy, grid_size=grid_size,
                                n_jobs=n_jobs,
                                basin_method=basin_method,
                                basin_connectivity=basin_connectivity,
                                basin_min_frames=basin_min_frames,
                                basin_min_depth=basin_min_depth,
                                mask_min_count=mask_min_count,
                                cmap=cmap, dpi=dpi)

        fel.main(energy_threshold, cv_names=cv_names, 
                 n_angles=n_angles, elevation=elevation, 
                 duration_per_frame=duration_per_frame)
        
        # Sempre grava a energia por frame. Antes so acontecia com --energy, mas a
        # tabela e o produto quantitativo principal e a coluna de bacia existe
        # independentemente de haver um corte de exibicao.
        print("Calculating and saving energy for each frame...")
        fel.calculate_and_save_free_energy(threshold=energy_threshold)
        print("Energy saved successfully!\n")

    except Exception as e:
        import traceback
        print(f"An unexpected error occurred: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
