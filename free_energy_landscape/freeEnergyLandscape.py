import os
import sys
import shutil
import platform
import tempfile
import subprocess
import numpy as np
import multiprocessing
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from joblib import Parallel, delayed
from matplotlib.colors import LinearSegmentedColormap

class FreeEnergyLandscape:

    def __init__(self, cv1_path, cv2_path, 
                 temperature, boltzmann_constant, 
                 bins=100, kde_bandwidth=None, 
                 cv_names=['CV1', 'CV2'], discrete=None,
                 xlim_inf=None, xlim_sup=None, ylim_inf=None, ylim_sup=None):  

        
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
        self.bins = bins
        self.kde_bandwidth = kde_bandwidth
        self.discrete = discrete
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
        # Carrega os dados das variáveis coletivas e os índices dos frames
        self.proj1_data_original = np.loadtxt(self.cv1_path, usecols=[1])
        self.proj2_data_original = np.loadtxt(self.cv2_path, usecols=[1])


    def calculate_free_energy(self, data):
        if hasattr(self, 'cached_results'):
            return self.cached_results

        values_original = np.vstack([data[:, 0], data[:, 1]]).T
        if self.kde_bandwidth:
            kernel_original = gaussian_kde(values_original.T, bw_method=self.kde_bandwidth)
        else:
            kernel_original = gaussian_kde(values_original.T)

        # Ajusta a geração da grade para respeitar os limites especificados
        x_min = self.xlim_inf if self.xlim_inf is not None else data[:, 0].min()
        x_max = self.xlim_sup if self.xlim_sup is not None else data[:, 0].max()
        y_min = self.ylim_inf if self.ylim_inf is not None else data[:, 1].min()
        y_max = self.ylim_sup if self.ylim_sup is not None else data[:, 1].max()
        
        X_original, Y_original = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]
        positions_original = np.vstack([X_original.ravel(), Y_original.ravel()])
        Z_original = np.reshape(kernel_original(positions_original).T, X_original.shape)
        G_original = -self.kB * self.temperature * np.log(Z_original)
        G_original = np.clip(G_original - np.min(G_original), 0, 25)

        self.cached_results = {'X_original': X_original, 
                            'Y_original': Y_original, 
                            'G_original': G_original
                            }
        return self.cached_results


    def boltzmann_inversion(self, data_list, titles, threshold=None):

        fig_combined, axs_combined = plt.subplots(1, len(data_list), figsize=(20, 6), sharey=True)        
        # Inicializa a lista de elementos da legenda com o elemento para a linha de energia livre
        legend_elements = [plt.Line2D([0], [0], color='red', lw=2, label='Free energy')]
        
        for idx, (ax, data, title) in enumerate(zip(axs_combined, data_list, titles)):
            hist, bin_edges = np.histogram(data, bins=self.bins, density=True)
            hist = np.clip(hist, a_min=1e-10, a_max=None)
            G = -self.kB * self.temperature * np.log(hist)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            G_min_normalized = G - np.min(G)
            
            # Plotar a linha de energia livre
            ax.plot(bin_centers, G_min_normalized, 
                    color='red', label='Free energy' 
                    if idx == 0 else "_nolegend_"
                    )
            
            ax.set_xlabel(title)
            
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
                    hist, bin_edges = np.histogram(data, bins=self.bins, density=True)
                    hist = np.clip(hist, a_min=1e-10, a_max=None)
                    G = -self.kB * self.temperature * np.log(hist)
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                    G_min_normalized = G - np.min(G)
                    mask = (G_min_normalized >= interval) & (G_min_normalized < end)
                    if np.any(mask):
                        ax.scatter(bin_centers[mask], G_min_normalized[mask], 
                                   color=self.discreet_colors[i % len(self.discreet_colors)], 
                                   marker=self.discreet_markers[i % len(self.discreet_markers)], 
                                   s=50)

        # Adiciona a legenda no último subplot
        axs_combined[-1].legend(handles=legend_elements, 
                                loc='upper left', 
                                bbox_to_anchor=(1.05, 1), 
                                title="Energy intervals")
        
        axs_combined[0].set_ylabel('Free Energy (kJ/mol)')
        plt.suptitle('Normalized Free Energy Profile Comparison')
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        plt.savefig('Combined_Free_Energy_Profile_Normalized.png')
        plt.show()
        

    def plot_histogram(self, data_list, titles):
        plt.figure(figsize=(8 * len(data_list), 6))
        
        # Normalizar os dados e calcular as frequências em porcentagem
        all_counts = []
        for data in data_list:
            data_normalized = (data - np.min(data)) / (np.max(data) - np.min(data))
            counts, _ = np.histogram(data_normalized, bins=self.bins)
            all_counts.append(counts)
        total_counts_max = max([max(counts) for counts in all_counts])
        
        for i, (data, title) in enumerate(zip(data_list, titles)):
            data_normalized = (data - np.min(data)) / (np.max(data) - np.min(data))
            counts, bin_edges = np.histogram(data_normalized, bins=self.bins)
            # Ajustar as contagens para que o valor máximo de contagem represente 100%
            normalized_counts = (counts / total_counts_max) * 100
            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            
            ax = plt.subplot(1, len(data_list), i + 1)
            bars = ax.bar(bin_centers, 
                          normalized_counts, 
                          width=(bin_centers[1] - bin_centers[0]), 
                          alpha=0.7, 
                          color='green'
                          )
            
            ax.set_ylim(0, 100)  # Definir explicitamente o eixo Y de 0 a 100%
            ax.set_title(f'Normalized {title}')
            ax.set_xlabel('Value')
            if i == 0:
                ax.set_ylabel('Frequency (%)')
            ax.grid(True)
        
        # Ajustar as marcas do eixo Y para incluir o valor máximo (100%)
        plt.yticks(np.arange(0, 101, 20))
        
        # Adicionar a legenda fora do loop, apenas uma vez
        plt.figlegend(['Normalized Frequency (%)'], loc='upper right')
        plt.tight_layout()
        plt.savefig('histograms_normalized_side_by_side.png')
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
        plt.savefig('cv_by_frame_combined_normalized.png')
        plt.show()


    def plot_energy_landscape(self, threshold, titles=['CV1', 'CV2'], xlim_inf=None, xlim_sup=None, ylim_inf=None, ylim_sup=None):

        if xlim_inf is not None and xlim_sup is not None:
            plt.xlim(xlim_inf, xlim_sup)
        if ylim_inf is not None and ylim_sup is not None:
            plt.ylim(ylim_inf, ylim_sup)

        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)
        plt.figure(figsize=(8, 6))

        custom_cmap = LinearSegmentedColormap.from_list("custom_energy", self.colors)

        # Define os níveis de contorno para uma variação suave de 2 em 2
        G_min, G_max = np.min(result['G_original']), np.max(result['G_original'])
        levels = np.arange(G_min, G_max, 2)
        
        cont = plt.contourf(result['X_original'], result['Y_original'], result['G_original'],
                            levels=levels, cmap=custom_cmap, extend='both')

        # Adiciona linhas de contorno para definição
        plt.contour(result['X_original'], result['Y_original'], result['G_original'],
                    levels=levels, colors='k', linewidths=0.5)

        # Lógica para plotar pontos discretizados se self.discrete for especificado
        if self.discrete is not None and threshold is not None:
            discrete_intervals = np.arange(0, threshold, self.discrete)

            for i, interval in enumerate(discrete_intervals):
                end = min(interval + self.discrete, threshold)
                mask = (result['G_original'].flatten() <= end) & (result['G_original'].flatten() > interval)
                X_flat, Y_flat = result['X_original'].flatten(), result['Y_original'].flatten()

                if np.any(mask):
                    plt.scatter(X_flat[mask], 
                                Y_flat[mask], 
                                color=self.discreet_colors[i % len(self.discreet_colors)],
                                marker=self.discreet_markers[i % len(self.discreet_markers)], 
                                label=f'{interval:.1f}-{end:.1f} KJ/mol'
                                )

        if threshold is not None:
            plt.legend(loc='lower left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

        cbar = plt.colorbar(cont)
        cbar.set_label('Free energy (kJ/mol)')
        plt.xlim(self.xlim_inf, self.xlim_sup)
        plt.ylim(self.ylim_inf, self.ylim_sup)
        plt.xlabel(titles[0])
        plt.ylabel(titles[1])
        plt.title('Free Energy Landscape')
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Ajuste para garantir que a legenda fique visível e não sobreponha o gráfico
        plt.savefig('Free_energy_landscape_with_discrete_points.png')
        plt.show()

    def plot_3D_energy_landscape(self, threshold=None, titles=['CV1', 'CV2'], xlim_inf=None, xlim_sup=None, ylim_inf=None, ylim_sup=None):
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)

        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')

        surf = ax.plot_surface(result['X_original'], 
                               result['Y_original'], 
                               result['G_original'], 
                               cmap=self.custom_cmap, 
                               edgecolor='none', 
                               alpha=0.6
                               )
        
        ax.set_xlabel(titles[0])
        ax.set_ylabel(titles[1])
        ax.set_zlabel('Free energy (kJ/mol)')
        ax.set_title('3D Free Energy Landscape')

        # Cores e marcadores para pontos discretizados
        if self.discrete is not None and threshold is not None:
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

        cbar = fig.colorbar(surf, 
                            shrink=0.5, 
                            aspect=5, 
                            pad=0.1, 
                            label='Free energy (kJ/mol)'
                            )
        
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), title="Energy intervals")
        
        plt.xlim(self.xlim_inf, self.xlim_sup)
        plt.ylim(self.ylim_inf, self.ylim_sup)
        plt.tight_layout()
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

        # Gera uma lista de ângulos para um movimento contínuo e suave
        angles = list(range(0, 360, int(360 / n_angles)))
        

        for i, angle in enumerate(angles):
            fig = plt.figure(figsize=(10, 7))
            ax = fig.add_subplot(111, projection='3d')

            surf = ax.plot_surface(
                result['X_original'], 
                result['Y_original'], 
                result['G_original'], 
                cmap=self.custom_cmap, 
                edgecolor='none', alpha=0.8, 
                vmin=np.min(result['G_original']), 
                vmax=np.max(result['G_original']))
            
            ax.view_init(elev=elevation, azim=angle)
            ax.set_xlabel(titles[0])
            ax.set_ylabel(titles[1])
            ax.set_zlabel('Free energy (kJ/mol)')
            ax.set_title('3D Free Energy Landscape')

            # Adiciona a barra de cores no primeiro e último frame
            if i == 0 or i == len(angles) - 1:
                fig.colorbar(surf, shrink=0.5, aspect=5, label='Free energy (kJ/mol)')

            frame_filename = os.path.join(temp_dir, f"frame_{angle:03d}.png")
            plt.savefig(frame_filename)
            filenames.append(frame_filename)
            plt.close()

        with imageio.get_writer(gif_filename, mode='I', duration=duration_per_frame) as writer:
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
            --energy                [float]     Energy (KJ/mol), single value (default: None)
            --discretize            [float]     Discrete value (KJ/mol) for energy (default: None)
            --bins_energy_histogram [int]       Bins for energy histogram (default: 100)
            --kde_bandwidth         [float]     Bandwidth for kernel density estimation (default: None)
            --names                 [str] [str] Names for the collective variables (default: CV1, CV2)
            --gif_angles            [int]       Angles for 3D GIF rotation (default: 10)
            --gif_elevation         [int]       Elevation angle for the 3D GIF (default: 10)
            --gif_duration          [float]     Duration per frame in the GIF in seconds (default: 0.1)
            --xlim_inf              [float]     Lower limit for the x-axis (default: None)
            --xlim_sup              [float]     Upper limit for the x-axis (default: None)
            --ylim_inf              [float]     Lower limit for the y-axis (default: None)
            --ylim_sup              [float]     Upper limit for the y-axis (default: None)

        Example:
            freeEnergyLandscape.py proj1Out.txt proj2Out.txt --energy 0.1 --discretize 0.1 --xlim_inf 27 --xlim_sup 130 --ylim_inf 3 --ylim_sup 37 

        """
        print(help_text)

    def calculate_density_for_chunk(self, combined_data_chunk, bw_method):
        # Esta função é uma versão simplificada que recalcula o kernel para cada chunk
        kernel = gaussian_kde(combined_data_chunk.T, bw_method=bw_method)
        density = np.exp(kernel.logpdf(combined_data_chunk.T))
        return density

    def calculate_and_save_free_energy(self, threshold=None):

        # Verifica se os dados foram carregados
        if self.proj1_data_original is None or self.proj2_data_original is None:
            raise ValueError("Data not loaded. Run load_data first.")

        # Carrega os índices dos frames
        frames = np.loadtxt(self.cv1_path, usecols=[0], dtype=np.float64).astype(np.int64)

        # Prepara os dados combinados
        combined_data = np.vstack((self.proj1_data_original, self.proj2_data_original)).T

        num_cpus = multiprocessing.cpu_count()
        data_chunks = np.array_split(combined_data, num_cpus, axis=0)

        # Recalcula a densidade de probabilidade para cada chunk de dados em paralelo
        results = Parallel(n_jobs=num_cpus)(delayed(self.calculate_density_for_chunk)(chunk, self.kde_bandwidth) for chunk in data_chunks)
        density = np.concatenate(results)
        
        # Calcula a energia livre
        G = -self.kB * self.temperature * np.log(density)
        G_min = np.min(G)
        G_normalized = G - G_min

        # Aplica o threshold, se especificado
        if threshold is not None:
            indices_below_threshold = G_normalized <= threshold
            filtered_frames = frames[indices_below_threshold]
            filtered_cv1 = self.proj1_data_original[indices_below_threshold]
            filtered_cv2 = self.proj2_data_original[indices_below_threshold]
            filtered_energy = G_normalized[indices_below_threshold]
        else:
            filtered_frames = frames
            filtered_cv1 = self.proj1_data_original
            filtered_cv2 = self.proj2_data_original
            filtered_energy = G_normalized

        # Prepara os dados para salvamento
        data_to_save = np.column_stack((filtered_frames, filtered_cv1, filtered_cv2, filtered_energy))

        # Ordena os dados pela energia
        data_to_save = data_to_save[data_to_save[:, 3].argsort()]

        # Salva os dados em um arquivo .tsv
        filename = 'discrete_values_energy_frames.tsv'
        np.savetxt(filename, 
                   data_to_save, 
                   delimiter='\t', 
                   fmt=['%d', '%.6f', '%.6f', '%.6f'], 
                   header='frame\tcv1\tcv2\tenergy', 
                   comments=''
                   )
        
        print(f"Energy data saved in'{filename}'.")


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

        print("Plotting the free energy landscape...\n")
        self.plot_energy_landscape(
            threshold=energy_threshold, titles=cv_names
            )
        print("Plot successfully generated.\n")

        print("Plotting the free energy landscape in 3D...\n")
        self.plot_3D_energy_landscape(
            threshold=energy_threshold, titles=cv_names
                                      )
        print("Plotting 3D gif...\n")

        self.create_3D_gif(
            n_angles=n_angles, elevation=elevation,
            duration_per_frame=duration_per_frame,
            titles=cv_names
                           )
        print("3D plot successfully generated.\n")

        # Após o uso final dos dados, limpe-os para liberar memória
        if hasattr(self, 'cached_results'):
            del self.cached_results

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
    xlim_inf = xlim_sup = ylim_inf = ylim_sup = None  # Inicialização padrão


    if len(sys.argv) >= 3:
        cv1_path, cv2_path = sys.argv[1], sys.argv[2]

        # Processar argumentos adicionais como pares chave-valor
        i = 3
        while i < len(sys.argv):
            key = sys.argv[i]
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
                                ylim_inf=ylim_inf, ylim_sup=ylim_sup)

        fel.main(energy_threshold, cv_names=cv_names, 
                 n_angles=n_angles, elevation=elevation, 
                 duration_per_frame=duration_per_frame)
        
        if energy_threshold is not None:
            print("Calculating and saving energy for each frame...")
            fel.calculate_and_save_free_energy(threshold=energy_threshold)
            print("Energy saved successfully!\n")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
