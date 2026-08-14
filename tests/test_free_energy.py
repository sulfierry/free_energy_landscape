"""Testes de corretude para o calculo de energia livre.

A ancora e analitica: uma mistura de duas gaussianas 2D bem separadas, com pesos
conhecidos w1 e w2, tem diferenca de energia livre entre os minimos igual a
    dG = -kB * T * ln(w1 / w2)
independente da largura das gaussianas (desde que estejam bem separadas).
"""

import numpy as np
import pytest
from scipy import ndimage

from free_energy_landscape.freeEnergyLandscape import FreeEnergyLandscape

KB = 8.314e-3   # kJ/(mol*K)
T = 300.0


def write_cv_files(tmp_path, cv1, cv2, frames=None):
    """Escreve dois arquivos no formato esperado: coluna 0 = frame, coluna 1 = valor."""
    if frames is None:
        frames = np.arange(len(cv1))
    p1 = tmp_path / "cv1.txt"
    p2 = tmp_path / "cv2.txt"
    np.savetxt(p1, np.column_stack([frames, cv1]))
    np.savetxt(p2, np.column_stack([frames, cv2]))
    return str(p1), str(p2)


def two_gaussian_mixture(w1, n=20000, seed=0, sigma=0.35, centers=((-3.0, 0.0), (3.0, 0.0))):
    """Amostra uma mistura de duas gaussianas 2D isotropicas com peso w1 na primeira."""
    rng = np.random.default_rng(seed)
    n1 = int(round(n * w1))
    n2 = n - n1
    a = rng.normal(centers[0], sigma, size=(n1, 2))
    b = rng.normal(centers[1], sigma, size=(n2, 2))
    data = np.vstack([a, b])
    rng.shuffle(data, axis=0)
    return data


# Duas gaussianas proximas o bastante para a regiao amostrada ficar conexa, de modo
# que exista uma sela real entre as bacias (em ±3 com sigma 0.35 elas ficam separadas
# por regiao sem amostragem nenhuma, e nao ha sela a medir).
ADJACENT = dict(sigma=0.6, centers=((-1.5, 0.0), (1.5, 0.0)))


def build(tmp_path, data, **kwargs):
    p1, p2 = write_cv_files(tmp_path, data[:, 0], data[:, 1])
    fel = FreeEnergyLandscape(p1, p2, T, KB, **kwargs)
    fel.load_data()
    return fel


# --------------------------------------------------------------------------
# Corretude fisica
# --------------------------------------------------------------------------

@pytest.mark.parametrize("w1, expected_dG", [
    (0.5, 0.0),
    (0.75, -KB * T * np.log(0.75 / 0.25)),
    (0.9, -KB * T * np.log(0.9 / 0.1)),
])
def test_recovers_analytical_delta_g(tmp_path, w1, expected_dG):
    """A FEL na grade reproduz o dG analitico entre os dois minimos."""
    data = two_gaussian_mixture(w1, n=30000, seed=1)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    res = fel.calculate_free_energy(data)

    X, Y, G = res['X_original'], res['Y_original'], res['G_original']

    # Minimo de G em cada bacia (x < 0 e x > 0). nanmin: fora da regiao
    # amostrada a superficie e NaN.
    left = np.nanmin(G[X < 0])
    right = np.nanmin(G[X > 0])
    measured_dG = left - right   # G(bacia 1) - G(bacia 2)

    assert measured_dG == pytest.approx(expected_dG, abs=0.35), (
        f"dG medido {measured_dG:.3f} kJ/mol, esperado {expected_dG:.3f} kJ/mol"
    )


def test_per_frame_energies_match_analytical_ratio(tmp_path):
    """A energia por frame (caminho do .tsv) tambem recupera o dG analitico."""
    w1 = 0.8
    data = two_gaussian_mixture(w1, n=20000, seed=2)
    fel = build(tmp_path, data, max_energy=None)

    kernel = fel.get_kernel(data)
    G = fel._to_free_energy(fel.evaluate_density(kernel, data))

    left = G[data[:, 0] < 0].min()
    right = G[data[:, 0] > 0].min()
    expected = -KB * T * np.log(w1 / (1 - w1))

    assert (left - right) == pytest.approx(expected, abs=0.35)


# --------------------------------------------------------------------------
# Regressao do bug critico: dependencia do numero de CPUs
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n_jobs", [1, 2, 3, 4, 8])
def test_density_is_independent_of_n_jobs(tmp_path, n_jobs):
    """A paralelizacao divide apenas a avaliacao: o resultado deve ser bit-identico.

    Regressao para o bug em que um gaussian_kde era reajustado dentro de cada chunk,
    tornando as densidades incomparaveis entre chunks e fazendo o resultado depender
    de multiprocessing.cpu_count().
    """
    data = two_gaussian_mixture(0.7, n=4000, seed=3)
    fel = build(tmp_path, data, max_energy=None)

    kernel = fel.get_kernel(data)
    serial = kernel(data.T)
    parallel = fel.evaluate_density(kernel, data, n_jobs=n_jobs)

    np.testing.assert_allclose(parallel, serial, rtol=1e-12, atol=0.0)


def test_selected_frames_are_independent_of_n_jobs(tmp_path):
    """O conjunto de frames abaixo do threshold nao pode mudar com o numero de CPUs."""
    data = two_gaussian_mixture(0.7, n=4000, seed=4)
    threshold = 2.0

    selections = []
    for n_jobs in (1, 2, 5, 8):
        fel = build(tmp_path, data, max_energy=None, n_jobs=n_jobs)
        out = tmp_path / f"energy_{n_jobs}.tsv"
        fel.calculate_and_save_free_energy(threshold=threshold, filename=str(out))
        selections.append(np.loadtxt(out, skiprows=1))

    for other in selections[1:]:
        np.testing.assert_allclose(other, selections[0], rtol=1e-10, atol=1e-10)


def test_chunked_refit_would_fail_this_test(tmp_path):
    """Documenta o comportamento antigo: reajustar o kernel por chunk diverge de verdade."""
    from scipy.stats import gaussian_kde

    data = two_gaussian_mixture(0.7, n=4000, seed=5)
    correct = gaussian_kde(data.T)(data.T)

    legacy = np.concatenate([
        np.exp(gaussian_kde(c.T).logpdf(c.T)) for c in np.array_split(data, 8, axis=0)
    ])

    # O metodo antigo NAO reproduz o kernel global - por isso foi removido.
    assert not np.allclose(legacy, correct, rtol=1e-3)


# --------------------------------------------------------------------------
# Cache
# --------------------------------------------------------------------------

def test_cache_is_invalidated_by_new_data(tmp_path):
    """calculate_free_energy nao pode devolver o resultado antigo para dados novos."""
    data_a = two_gaussian_mixture(0.5, n=3000, seed=6)
    data_b = two_gaussian_mixture(0.5, n=3000, seed=7,
                                  centers=((-8.0, -8.0), (8.0, 8.0)))

    fel = build(tmp_path, data_a, max_energy=None)
    first = fel.calculate_free_energy(data_a)['X_original'].copy()
    second = fel.calculate_free_energy(data_b)['X_original']

    assert not np.allclose(first, second), "cache devolveu a grade antiga para dados novos"


def test_cache_hit_returns_same_object(tmp_path):
    data = two_gaussian_mixture(0.5, n=3000, seed=8)
    fel = build(tmp_path, data, max_energy=None)
    assert fel.calculate_free_energy(data) is fel.calculate_free_energy(data)


def test_cache_is_invalidated_by_bandwidth(tmp_path):
    data = two_gaussian_mixture(0.5, n=3000, seed=9)
    fel = build(tmp_path, data, max_energy=None)
    coarse = fel.calculate_free_energy(data)['G_original'].copy()
    fel.kde_bandwidth = 0.15
    fine = fel.calculate_free_energy(data)['G_original']
    assert not np.allclose(coarse, fine, equal_nan=True)


# --------------------------------------------------------------------------
# max_energy
# --------------------------------------------------------------------------

def test_max_energy_none_preserves_the_full_range(tmp_path):
    """Sem cap, a barreira real sobrevive; com cap, e truncada."""
    data = two_gaussian_mixture(0.5, n=20000, seed=10, sigma=0.3)
    uncapped = build(tmp_path, data, max_energy=None, grid_size=100)
    full = np.nanmax(uncapped.calculate_free_energy(data)['G_original'])

    cap = full / 2
    capped = build(tmp_path, data, max_energy=cap, grid_size=100)
    G = capped.calculate_free_energy(data)['G_original']

    assert full > cap
    assert np.nanmax(G) == pytest.approx(cap)
    assert np.nanmin(G) == pytest.approx(0.0)


def test_unsampled_region_is_masked_out(tmp_path):
    """O KDE nao pode pintar energia onde nao ha um unico frame."""
    data = two_gaussian_mixture(0.5, n=20000, seed=10, sigma=0.3)
    fel = build(tmp_path, data, max_energy=None, grid_size=100)
    res = fel.calculate_free_energy(data)

    assert res['mask'].any() and not res['mask'].all()
    assert np.all(np.isnan(res['G_original'][~res['mask']]))
    assert np.all(np.isfinite(res['G_original'][res['mask']]))
    # o interior da mascara e onde o KDE espera pelo menos um frame por celula
    interior = ndimage.binary_erosion(res['mask'], np.ones((3, 3), bool))
    assert res['expected_counts'][interior].min() >= 1.0


def test_mask_min_count_controls_the_extent(tmp_path):
    data = two_gaussian_mixture(0.5, n=20000, seed=10, sigma=0.3)
    loose = build(tmp_path, data, max_energy=None, mask_min_count=0.1)
    strict = build(tmp_path, data, max_energy=None, mask_min_count=10.0)
    assert (loose.calculate_free_energy(data)['mask'].sum()
            > strict.calculate_free_energy(data)['mask'].sum())


def test_cap_does_not_leak_into_per_frame_energies(tmp_path):
    """O cap e so escala de cor: o .tsv nunca deve ser truncado."""
    data = two_gaussian_mixture(0.5, n=8000, seed=11, sigma=0.3)
    fel = build(tmp_path, data, max_energy=1.0)
    out = tmp_path / "energy.tsv"
    fel.calculate_and_save_free_energy(filename=str(out))
    saved = np.loadtxt(out, skiprows=1)
    assert saved[:, 3].max() > 1.0


# --------------------------------------------------------------------------
# Validacao de entrada
# --------------------------------------------------------------------------

def test_rejects_mismatched_lengths(tmp_path):
    p1 = tmp_path / "a.txt"
    p2 = tmp_path / "b.txt"
    np.savetxt(p1, np.column_stack([np.arange(10), np.zeros(10)]))
    np.savetxt(p2, np.column_stack([np.arange(8), np.zeros(8)]))
    fel = FreeEnergyLandscape(str(p1), str(p2), T, KB)
    with pytest.raises(ValueError, match="different numbers of frames"):
        fel.load_data()


def test_rejects_mismatched_frame_indices(tmp_path):
    p1 = tmp_path / "a.txt"
    p2 = tmp_path / "b.txt"
    np.savetxt(p1, np.column_stack([np.arange(10), np.zeros(10)]))
    np.savetxt(p2, np.column_stack([np.arange(10) + 5, np.zeros(10)]))
    fel = FreeEnergyLandscape(str(p1), str(p2), T, KB)
    with pytest.raises(ValueError, match="Frame indices do not match"):
        fel.load_data()


def test_rejects_single_column_file(tmp_path):
    p1 = tmp_path / "a.txt"
    p2 = tmp_path / "b.txt"
    np.savetxt(p1, np.zeros(10))
    np.savetxt(p2, np.column_stack([np.arange(10), np.zeros(10)]))
    fel = FreeEnergyLandscape(str(p1), str(p2), T, KB)
    with pytest.raises(ValueError, match="at least 2 columns"):
        fel.load_data()


def test_rejects_non_finite_values(tmp_path):
    cv = np.zeros(10)
    cv[3] = np.nan
    p1, p2 = write_cv_files(tmp_path, cv, np.zeros(10))
    fel = FreeEnergyLandscape(p1, p2, T, KB)
    with pytest.raises(ValueError, match="NaN or infinite"):
        fel.load_data()


def test_load_data_preserves_frame_indices(tmp_path):
    frames = np.arange(100, 110)
    data = two_gaussian_mixture(0.5, n=10, seed=12)
    p1, p2 = write_cv_files(tmp_path, data[:, 0], data[:, 1], frames=frames)
    fel = FreeEnergyLandscape(p1, p2, T, KB)
    fel.load_data()
    np.testing.assert_array_equal(fel.frames, frames)


def test_saved_tsv_uses_original_frame_indices(tmp_path):
    frames = np.arange(1000, 1000 + 2000)
    data = two_gaussian_mixture(0.6, n=2000, seed=13)
    p1, p2 = write_cv_files(tmp_path, data[:, 0], data[:, 1], frames=frames)
    fel = FreeEnergyLandscape(p1, p2, T, KB, max_energy=None)
    fel.load_data()
    out = tmp_path / "energy.tsv"
    fel.calculate_and_save_free_energy(filename=str(out))
    saved = np.loadtxt(out, skiprows=1)
    assert saved.shape[0] == len(frames)
    assert set(saved[:, 0].astype(int)) == set(frames.tolist())


def test_tsv_is_sorted_by_energy(tmp_path):
    data = two_gaussian_mixture(0.6, n=2000, seed=14)
    fel = build(tmp_path, data, max_energy=None)
    out = tmp_path / "energy.tsv"
    fel.calculate_and_save_free_energy(filename=str(out))
    energy = np.loadtxt(out, skiprows=1)[:, 3]
    assert np.all(np.diff(energy) >= 0)


def test_threshold_filters_frames(tmp_path):
    data = two_gaussian_mixture(0.6, n=3000, seed=15)
    fel = build(tmp_path, data, max_energy=None)
    full = tmp_path / "full.tsv"
    cut = tmp_path / "cut.tsv"
    fel.calculate_and_save_free_energy(filename=str(full))
    fel.calculate_and_save_free_energy(threshold=1.0, filename=str(cut))
    a = np.loadtxt(full, skiprows=1)
    b = np.loadtxt(cut, skiprows=1)
    assert b.shape[0] < a.shape[0]
    assert b[:, 3].max() <= 1.0


def test_requires_load_data_first(tmp_path):
    data = two_gaussian_mixture(0.5, n=100, seed=16)
    p1, p2 = write_cv_files(tmp_path, data[:, 0], data[:, 1])
    fel = FreeEnergyLandscape(p1, p2, T, KB)
    with pytest.raises(ValueError, match="Data not loaded"):
        fel.calculate_and_save_free_energy()


# --------------------------------------------------------------------------
# Clusterizacao em bacias
# --------------------------------------------------------------------------

def test_finds_exactly_two_basins(tmp_path):
    """Duas gaussianas bem separadas devem dar duas bacias, nao dois grupos de pontos soltos."""
    data = two_gaussian_mixture(0.7, n=20000, seed=20)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    result = fel.identify_basins(threshold=3.0)
    assert len(result['basins']) == 2


def test_basins_are_ordered_by_energy(tmp_path):
    data = two_gaussian_mixture(0.7, n=20000, seed=21)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    basins = fel.identify_basins(threshold=3.0)['basins']
    assert [b['id'] for b in basins] == list(range(1, len(basins) + 1))
    assert basins[0]['G_min'] == pytest.approx(0.0, abs=1e-9)
    assert all(a['G_min'] <= b['G_min'] for a, b in zip(basins, basins[1:]))


def test_basin_minima_reproduce_analytical_delta_g(tmp_path):
    """A diferenca de energia entre as bacias e a resposta analitica da mistura."""
    w1 = 0.75
    data = two_gaussian_mixture(w1, n=30000, seed=22)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    basins = fel.identify_basins(threshold=4.0)['basins']
    expected = -KB * T * np.log(w1 / (1 - w1))
    assert basins[1]['G_min'] - basins[0]['G_min'] == pytest.approx(abs(expected), abs=0.35)


def test_basin_minima_sit_on_the_true_centres(tmp_path):
    data = two_gaussian_mixture(0.6, n=20000, seed=23)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    basins = fel.identify_basins(threshold=3.0)['basins']
    centres = sorted(b['cv1_min'] for b in basins)
    assert centres[0] == pytest.approx(-3.0, abs=0.4)
    assert centres[1] == pytest.approx(3.0, abs=0.4)
    for b in basins:
        assert b['cv2_min'] == pytest.approx(0.0, abs=0.4)


def test_representative_frame_is_a_real_frame(tmp_path):
    """O ponto representativo e uma conformacao amostrada, nao um no da grade."""
    frames = np.arange(500, 500 + 20000)
    data = two_gaussian_mixture(0.6, n=20000, seed=24)
    p1, p2 = write_cv_files(tmp_path, data[:, 0], data[:, 1], frames=frames)
    fel = FreeEnergyLandscape(p1, p2, T, KB, max_energy=None, grid_size=120)
    fel.load_data()
    basins = fel.identify_basins(threshold=3.0)['basins']

    for b in basins:
        assert b['rep_frame'] in set(frames.tolist())
        row = int(np.flatnonzero(frames == b['rep_frame'])[0])
        assert b['rep_cv1'] == pytest.approx(data[row, 0])
        assert b['rep_cv2'] == pytest.approx(data[row, 1])
        # o representante esta perto do minimo da sua bacia
        assert abs(b['rep_cv1'] - b['cv1_min']) < 1.0


def test_populations_sum_to_the_frames_below_threshold(tmp_path):
    data = two_gaussian_mixture(0.6, n=10000, seed=25)
    fel = build(tmp_path, data, max_energy=None, grid_size=100)
    result = fel.identify_basins(threshold=3.0)
    assigned = sum(b['n_frames'] for b in result['basins'])
    assert assigned == int(np.count_nonzero(result['frame_basin'] > 0))
    assert 0 < assigned <= len(data)


def test_frame_basin_matches_grid_labels(tmp_path):
    data = two_gaussian_mixture(0.6, n=5000, seed=26)
    fel = build(tmp_path, data, max_energy=None, grid_size=100)
    result = fel.identify_basins(threshold=3.0)
    grid = fel.calculate_free_energy(data)
    i, j, inside = fel._grid_indices(grid)
    expected = np.where(inside, result['labels'][i, j], 0)
    np.testing.assert_array_equal(result['frame_basin'], expected)


def test_min_frames_discards_small_basins(tmp_path):
    data = two_gaussian_mixture(0.7, n=20000, seed=27)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    all_basins = fel.identify_basins(threshold=3.0)['basins']
    biggest = max(b['n_frames'] for b in all_basins)
    filtered = fel.identify_basins(threshold=3.0, min_frames=biggest)['basins']
    assert len(filtered) == 1
    assert len(filtered) < len(all_basins)


def test_threshold_above_cap_is_rejected(tmp_path):
    """Acima do cap a superficie e constante e todas as bacias se fundiriam numa so."""
    data = two_gaussian_mixture(0.5, n=5000, seed=28)
    fel = build(tmp_path, data, max_energy=10.0)
    with pytest.raises(ValueError, match="above --max_energy"):
        fel.identify_basins(threshold=20.0, method='connected')


def test_connected_merges_basins_at_higher_threshold(tmp_path):
    """No modo 'connected', subir o corte acima da barreira funde as duas bacias."""
    data = two_gaussian_mixture(0.5, n=30000, seed=29, **ADJACENT)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    saddle = min(fel._basin_saddles(
        fel.calculate_free_energy(data)['G_original'],
        fel.identify_basins(min_depth=0.0)['labels']).values())
    assert len(fel.identify_basins(threshold=saddle * 0.9,
                                   method='connected')['basins']) == 2
    assert len(fel.identify_basins(threshold=saddle * 1.1,
                                   method='connected')['basins']) == 1


def test_connected_requires_a_threshold(tmp_path):
    data = two_gaussian_mixture(0.5, n=3000, seed=33)
    fel = build(tmp_path, data, max_energy=None)
    with pytest.raises(ValueError, match="requires an energy threshold"):
        fel.identify_basins(method='connected')


def test_unknown_method_is_rejected(tmp_path):
    data = two_gaussian_mixture(0.5, n=3000, seed=34)
    fel = build(tmp_path, data, max_energy=None)
    with pytest.raises(ValueError, match="Unknown basin method"):
        fel.identify_basins(method='kmeans')


# --------------------------------------------------------------------------
# Watershed: bacia de atracao
# --------------------------------------------------------------------------

def test_watershed_assigns_every_sampled_frame(tmp_path):
    """Ao contrario de 'connected', o watershed nao deixa frame orfao na regiao amostrada."""
    data = two_gaussian_mixture(0.6, n=20000, seed=35)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    result = fel.identify_basins()

    res = fel.calculate_free_energy(data)
    i, j, inside = fel._grid_indices(res)
    sampled = inside & res['mask'][i, j]
    assert np.all(result['frame_basin'][sampled] > 0)


def test_watershed_needs_no_threshold(tmp_path):
    data = two_gaussian_mixture(0.6, n=20000, seed=36)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    assert len(fel.identify_basins()['basins']) >= 1


def test_watershed_threshold_only_filters_reporting(tmp_path):
    """O threshold no watershed nao muda o particionamento, so o que e reportado."""
    data = two_gaussian_mixture(0.75, n=30000, seed=37)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    everything = fel.identify_basins()['basins']
    assert len(everything) >= 2

    cut = everything[0]['G_min'] + 0.5 * (everything[1]['G_min'] - everything[0]['G_min'])
    filtered = fel.identify_basins(threshold=cut)['basins']
    assert len(filtered) == 1
    assert filtered[0]['n_frames'] == everything[0]['n_frames']


def test_basin_depth_is_the_saddle_height(tmp_path):
    """A profundidade e a altura da sela mais baixa acima do minimo da bacia."""
    data = two_gaussian_mixture(0.5, n=30000, seed=38, **ADJACENT)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    result = fel.identify_basins(min_depth=0.0)
    basins = result['basins']
    assert len(basins) == 2

    G = fel.calculate_free_energy(data)['G_original']
    saddles = fel._basin_saddles(G, result['labels'])
    assert saddles, "as duas bacias deveriam ser adjacentes"
    height = min(saddles.values())
    for b in basins:
        assert b['depth'] == pytest.approx(height - b['G_min'], abs=1e-9)


def test_min_depth_merges_shallow_basins(tmp_path):
    """Barreiras rasas sao ruido do KDE: min_depth alto colapsa tudo numa bacia."""
    data = two_gaussian_mixture(0.5, n=30000, seed=39, **ADJACENT)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    fine = fel.identify_basins(min_depth=0.0)['basins']
    assert len(fine) == 2

    barrier = max(b['depth'] for b in fine)
    coarse = fel.identify_basins(min_depth=barrier + 1.0)['basins']
    assert len(coarse) == 1
    assert coarse[0]['n_frames'] > fine[0]['n_frames']


def test_default_min_depth_is_auto(tmp_path):
    data = two_gaussian_mixture(0.5, n=3000, seed=40)
    fel = build(tmp_path, data, max_energy=None)
    assert fel.basin_min_depth == 'auto'


def test_auto_min_depth_cuts_at_the_largest_gap(tmp_path):
    """O corte automatico cai no maior salto do espectro de persistencia."""
    data = two_gaussian_mixture(0.5, n=3000, seed=42)
    fel = build(tmp_path, data, max_energy=None)

    # espectro sintetico: estrutura real bem acima do ruido
    assert fel._auto_min_depth([2.1, 1.8, 1.4, 0.6, 0.4, 0.1]) == pytest.approx(1.0)
    # nunca acima de kB*T, mesmo com um salto enorme no topo
    assert fel._auto_min_depth([40.0, 0.1]) == pytest.approx(KB * T)
    # sem espectro nao ha o que fundir
    assert fel._auto_min_depth([]) == 0.0
    assert fel._auto_min_depth([1.0]) == 0.0


def test_auto_min_depth_is_reported(tmp_path):
    data = two_gaussian_mixture(0.6, n=20000, seed=43)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    result = fel.identify_basins()
    assert result['auto_depth'] is True
    assert result['n_raw_minima'] >= len(result['basins'])
    assert isinstance(result['min_depth'], float)


def test_explicit_min_depth_overrides_auto(tmp_path):
    data = two_gaussian_mixture(0.6, n=20000, seed=44)
    fel = build(tmp_path, data, max_energy=None, grid_size=120, basin_min_depth=0.0)
    result = fel.identify_basins()
    assert result['auto_depth'] is False
    assert result['min_depth'] == 0.0


def test_watershed_populations_cover_the_sampled_frames(tmp_path):
    data = two_gaussian_mixture(0.7, n=20000, seed=41)
    fel = build(tmp_path, data, max_energy=None, grid_size=120)
    result = fel.identify_basins()
    assigned = sum(b['n_frames'] for b in result['basins'])
    assert assigned == int(np.count_nonzero(result['frame_basin'] > 0))
    # a esmagadora maioria dos frames cai dentro da regiao amostrada
    assert assigned > 0.95 * len(data)


def test_tsv_gains_a_cluster_column(tmp_path):
    data = two_gaussian_mixture(0.6, n=5000, seed=30)
    fel = build(tmp_path, data, max_energy=None, grid_size=100)
    fel._basins = fel.identify_basins(threshold=3.0)
    out = tmp_path / "energy.tsv"
    fel.calculate_and_save_free_energy(threshold=3.0, filename=str(out))

    with open(out) as fh:
        header = fh.readline().strip().split('\t')
    assert header == ['frame', 'cv1', 'cv2', 'energy', 'cluster']

    saved = np.loadtxt(out, skiprows=1)
    assert set(np.unique(saved[:, 4]).astype(int)) <= {0, 1, 2}


def test_save_basins_roundtrip(tmp_path):
    data = two_gaussian_mixture(0.6, n=5000, seed=31)
    fel = build(tmp_path, data, max_energy=None, grid_size=100)
    basins = fel.identify_basins(threshold=3.0)['basins']
    out = tmp_path / "basins.tsv"
    fel.save_basins(basins, filename=str(out))

    with open(out) as fh:
        header = fh.readline().strip().split('\t')
    assert header[0] == 'cluster' and 'rep_frame' in header

    saved = np.atleast_2d(np.loadtxt(out, skiprows=1))
    assert saved.shape[0] == len(basins)
    assert saved[0, 0] == 1


def test_basins_are_independent_of_n_jobs(tmp_path):
    data = two_gaussian_mixture(0.7, n=6000, seed=32)
    reference = None
    for n_jobs in (1, 3, 8):
        fel = build(tmp_path, data, max_energy=None, grid_size=100, n_jobs=n_jobs)
        basins = fel.identify_basins(threshold=3.0)['basins']
        summary = [(b['id'], b['n_frames'], b['rep_frame'], round(b['G_min'], 9))
                   for b in basins]
        if reference is None:
            reference = summary
        assert summary == reference


# --------------------------------------------------------------------------
# Grade
# --------------------------------------------------------------------------

def test_grid_size_is_honoured(tmp_path):
    data = two_gaussian_mixture(0.5, n=2000, seed=17)
    fel = build(tmp_path, data, grid_size=37, max_energy=None)
    assert fel.calculate_free_energy(data)['G_original'].shape == (37, 37)


def test_axis_limits_define_the_grid(tmp_path):
    data = two_gaussian_mixture(0.5, n=20000, seed=18)
    fel = build(tmp_path, data, max_energy=None,
                xlim_inf=-5.0, xlim_sup=5.0, ylim_inf=-2.5, ylim_sup=2.5)
    res = fel.calculate_free_energy(data)
    assert res['X_original'].min() == pytest.approx(-5.0)
    assert res['X_original'].max() == pytest.approx(5.0)
    assert res['Y_original'].min() == pytest.approx(-2.5)
    assert res['Y_original'].max() == pytest.approx(2.5)


def test_view_limits_zoom_to_the_sampled_region(tmp_path):
    """Os eixos exibidos acompanham a regiao amostrada, nao o retangulo inteiro."""
    data = two_gaussian_mixture(0.5, n=20000, seed=19)
    fel = build(tmp_path, data, max_energy=None,
                xlim_inf=-20.0, xlim_sup=20.0, ylim_inf=-20.0, ylim_sup=20.0)
    res = fel.calculate_free_energy(data)
    (x_lo, x_hi), (y_lo, y_hi) = fel._view_limits(res)
    # os limites explicitos tem precedencia
    assert (x_lo, x_hi, y_lo, y_hi) == (-20.0, 20.0, -20.0, 20.0)

    fel.xlim_inf = fel.xlim_sup = fel.ylim_inf = fel.ylim_sup = None
    (x_lo, x_hi), (y_lo, y_hi) = fel._view_limits(res)
    assert -20.0 < x_lo < x_hi < 20.0
    assert -20.0 < y_lo < y_hi < 20.0
