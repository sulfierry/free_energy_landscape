import os
import sys
import shutil
import platform
import tempfile
import subprocess
import numpy as np
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap

class FreeEnergyLandscape:
    
    def __init__(self, cv1_path, cv2_path, temperature, boltzmann_constant, bins=100, kde_bandwidth=None, cv_names=['CV1', 'CV2']):
        self.cv1_path = cv1_path
        self.cv2_path = cv2_path
        self.temperature = temperature
        self.kB = boltzmann_constant
        self.cv_names = cv_names
        self.colors = [
            (0, "darkblue"),    # 0 a 3
            (3/25, "blue"),     # 3 a 6
            (6/25, "lightblue"),# 6 a 9
            (9/25, "#ADD8E6"),  # 9 a 12 azul claríssimo
            (12/25, "#FFA07A"), # 12 a 15 vermelho claro (quase salmão)
            (15/25, "#FF4500"), # 15 a 18 mais escuro (quase laranja)
            (18/25, "#FF6347"), # 18 a 21 laranja/vermelho
            (21/25, "darkred"), # 21 a 24 vermelho escuro
            (1, "darkred")      # Garante que o máximo seja vermelho escuro
        ]
        self.custom_cmap = LinearSegmentedColormap.from_list("custom_energy", self.colors)
        self.proj1_data_original = None
        self.proj2_data_original = None
        self.bins = bins
        self.kde_bandwidth = kde_bandwidth
        

    def load_data(self):
        self.proj1_data_original = np.loadtxt(self.cv1_path, usecols=[1])
        self.proj2_data_original = np.loadtxt(self.cv2_path, usecols=[1])


    def boltzmann_inversion(self, data_list, titles, threshold):
        fig_combined, axs_combined = plt.subplots(1, len(data_list), figsize=(20, 6), sharey=True)

        # Primeiro, calculamos o G_max_global com base em todos os conjuntos de dados
        G_max_global = -np.inf
        for data in data_list:
            hist, bin_edges = np.histogram(data, bins=100, density=True)
            hist = np.clip(hist, a_min=1e-10, a_max=None)
            G = -self.kB * self.temperature * np.log(hist + 1e-10)
            G_max_global = max(G_max_global, np.max(G))

        # Agora, plotamos cada variável coletiva e seus pontos de baixa energia
        for ax, data, title in zip(axs_combined, data_list, titles):
            hist, bin_edges = np.histogram(data, bins=100, density=True)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            G = -self.kB * self.temperature * np.log(hist + 1e-10)
            G_normalized = G / G_max_global

            ax.plot(bin_centers, G_normalized, label='Free energy', color='red')
            ax.set_xlabel(title)
            ax.set_ylim(0, 1)

            # Aplicamos o threshold diretamente ao G normalizado para cada variável
            if threshold is not None:
                # Para cada variável, recalculamos o threshold_normalized baseado no seu próprio G_max
                G_max = np.max(G)
                threshold_normalized = threshold / G_max
                low_energy_bins = G_normalized < threshold_normalized
                if any(low_energy_bins):
                    ax.scatter(bin_centers[low_energy_bins], G_normalized[low_energy_bins], color='magenta', label='Low energy points')

        axs_combined[0].set_ylabel('Normalized Free Energy')
        plt.legend()
        plt.suptitle('Normalized Free Energy Profile Comparison')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig('Combined_Free_Energy_Profile_Normalized.png')
        plt.show()

    def calculate_free_energy(self, data):
        if hasattr(self, 'cached_results'):
            return self.cached_results

        values_original = np.vstack([data[:, 0], data[:, 1]]).T
        if self.kde_bandwidth:
            kernel_original = gaussian_kde(values_original.T, bw_method=self.kde_bandwidth)
        else:
            kernel_original = gaussian_kde(values_original.T)
        X_original, Y_original = np.mgrid[data[:, 0].min():data[:, 0].max():100j, 
                                          data[:, 1].min():data[:, 1].max():100j]
        positions_original = np.vstack([X_original.ravel(), Y_original.ravel()])
        Z_original = np.reshape(kernel_original(positions_original).T, X_original.shape)
        G_original = -self.kB * self.temperature * np.log(Z_original)
        G_original = np.clip(G_original - np.min(G_original), 0, 25)

        self.cached_results = {'X_original': X_original, 'Y_original': Y_original, 'G_original': G_original}
        return self.cached_results


    def plot_energy_landscape(self, threshold, titles=['CV1', 'CV2']):
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)
        plt.figure(figsize=(8, 6))

        cont = plt.contourf(result['X_original'], result['Y_original'], result['G_original'], 
                            levels=np.linspace(np.min(result['G_original']), np.max(result['G_original']), 100), 
                            cmap=self.custom_cmap, extend='both')

        if threshold is not None:
            low_energy_mask = result['G_original'] <= threshold
            if np.any(low_energy_mask):  # Verifica se existem pontos de baixa energia para plotar
                plt.scatter(result['X_original'][low_energy_mask], result['Y_original'][low_energy_mask], 
                            color='magenta', s=10, label=f'Energy <= {threshold} kJ/mol')
                plt.legend(loc='lower left', bbox_to_anchor=(1, 1))

        if threshold is not None:
            with open('landscape_cv1_cv2_minimuns.tsv', 'w') as file:
                file.write("frame\tcv1\tcv2\tenergy\n")
                for i in np.where(low_energy_mask.ravel())[0]:
                    energy = result['G_original'].ravel()[i]
                    cv1 = result['X_original'].ravel()[i]
                    cv2 = result['Y_original'].ravel()[i]
                    file.write(f"{i}\t{cv1}\t{cv2}\t{energy}\n")

        cbar = plt.colorbar(cont)
        cbar.set_label('Free energy (kJ/mol)')
        # Usa os títulos fornecidos para os eixos
        plt.xlabel(titles[0])
        plt.ylabel(titles[1])
        plt.title('Free energy landscape')
        plt.savefig('Free_energy_landscape.png')
        plt.show()



    def save_low_energy_points_to_tsv(self, threshold):
        if threshold is None or threshold < 0:
            print("Threshold is None or less than 0. No points will be saved.")
            return

        # Carrega os dados originais novamente para obter os números dos frames
        cv1_data = np.loadtxt(self.cv1_path)
        cv2_data = np.loadtxt(self.cv2_path)

        # Supondo que a primeira coluna de cada arquivo contém os números dos frames
        frames_cv1 = cv1_data[:, 0].astype(int)
        frames_cv2 = cv2_data[:, 0].astype(int)

        # Calcula a paisagem de energia usando os dados carregados previamente
        result = self.calculate_free_energy(np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None])))

        low_energy_mask = result['G_original'] <= threshold
        if not np.any(low_energy_mask):
            print("No low energy points found. No points will be saved.")
            return

        with open('low_energy_points.tsv', 'w') as file:
            file.write("frame\tcv1\tcv2\tenergy\n")

            # Itera sobre os pontos de baixa energia
            for idx in np.where(low_energy_mask.ravel())[0]:
                # A coordenada X, Y no espaço de CV1 e CV2 pode ser mapeada de volta para o frame mais próximo
                # Aqui, assumimos que os frames em CV1 e CV2 estão alinhados e usamos frames de CV1 como referência
                frame = frames_cv1[idx]  # Usa o número do frame diretamente dos dados carregados
                cv1 = self.proj1_data_original[idx]
                cv2 = self.proj2_data_original[idx]
                g = result['G_original'].flatten()[idx]
                file.write(f"{frame}\t{cv1}\t{cv2}\t{g}\n")

    def plot_3D_energy_landscape(self, threshold=None, titles=['CV1', 'CV2']):
        data = np.hstack((self.proj1_data_original[:, None], self.proj2_data_original[:, None]))
        result = self.calculate_free_energy(data)
        
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plotar a superfície da paisagem de energia
        surf = ax.plot_surface(result['X_original'], result['Y_original'], result['G_original'], cmap=self.custom_cmap, edgecolor='none', alpha=0.5)
        ax.set_xlabel(titles[0])
        ax.set_ylabel(titles[1])
        ax.set_zlabel('Free energy (kJ/mol)')
        ax.set_title('3D Free Energy Landscape')
        
        if threshold is not None:
            if isinstance(threshold, list):
                for i, interval in enumerate(threshold):
                    self.plot_threshold_points(ax, result, interval[0], interval[1], self.colors[i % len(self.colors)], f'Energy {interval[0]}-{interval[1]} kJ/mol')
            elif isinstance(threshold, (int, float)):
                self.plot_threshold_points(ax, result, 0, threshold, 'magenta', f'Energy <= {threshold}')
        
        fig.colorbar(surf, shrink=0.5, aspect=5, label='Free energy (kJ/mol)')
        plt.legend(loc='upper left', frameon=False)  # Adicionado frameon=False para evitar erros caso não haja legendas a serem exibidas
        plt.show()


    def plot_threshold_points(self, ax, result, lower_bound, upper_bound, color, label):
        G_flat = result['G_original'].flatten()
        energy_mask = (G_flat >= lower_bound) & (G_flat < upper_bound)
        
        if any(energy_mask):
            X_flat, Y_flat = result['X_original'].flatten(), result['Y_original'].flatten()
            ax.scatter(X_flat[energy_mask], Y_flat[energy_mask], G_flat[energy_mask], color=color, s=20, label=label)

    def create_3D_gif(self, gif_filename='energy_landscape_3D.gif', n_angles=10, elevation=15, duration_per_frame=0.01, titles=['CV1', 'CV2']):
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
            surf = ax.plot_surface(result['X_original'], result['Y_original'], result['G_original'], cmap=self.custom_cmap, edgecolor='none', vmin=np.min(result['G_original']), vmax=np.max(result['G_original']))
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

    def plot_histogram(self, data_list, titles):
        # Plotando histogramas absolutos individualmente
        for data, title in zip(data_list, titles):
            plt.figure(figsize=(8, 6))
            plt.hist(data, bins=self.bins, density=True, alpha=0.7, color='blue')
            plt.title(f'Absolute {title}')
            plt.xlabel('Value')
            plt.ylabel('Frequency')
            plt.grid(True)
            # plt.savefig(f'histogram_absolute_{title.replace(" ", "_")}.png')
            plt.close()  # Fecha a figura após salvar

        # Plotando histogramas normalizados lado a lado
        plt.figure(figsize=(8 * len(data_list), 6))
        for i, (data, title) in enumerate(zip(data_list, titles)):
            # Normalização dos dados
            data_normalized = (data - np.min(data)) / (np.max(data) - np.min(data))
            
            # Criação de subplots lado a lado
            ax = plt.subplot(1, len(data_list), i + 1)
            ax.hist(data_normalized, bins=self.bins, density=True, alpha=0.7, color='green')
            ax.set_title(f'Normalized {title}')
            ax.set_xlabel('Normalized Value')
            if i == 0:  # Apenas o primeiro subplot terá o label do eixo Y
                ax.set_ylabel('Frequency')
            ax.grid(True)

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

        # Plot combinado dos valores relativos
        plt.figure(figsize=(10, 6))
        for data, title in zip(data_list, titles):
            data_normalized = (data - np.min(data)) / (np.max(data) - np.min(data))
            plt.plot(frames, data_normalized, label=title)
        plt.xlabel('Frame')
        plt.ylabel('Normalized CV')
        plt.title('CV by Frame - Combined Normalized')
        plt.legend()
        plt.savefig('cv_by_frame_combined_normalized.png')
        plt.show()


    def verify_input(self, data_path):
        try:
            data = np.loadtxt(data_path, dtype=float)
            if data.shape[1] != 2:
                raise ValueError("O arquivo de entrada deve ter exatamente duas colunas.")

            frames = data[:, 0]
            if not np.all(frames == np.arange(len(frames))):
                raise ValueError("A primeira coluna deve conter valores inteiros sequenciais que representam os frames.")

        except ValueError as e:
            print(f"Erro na verificação do arquivo {data_path}: {e}")
            sys.exit(1)

    @staticmethod
    def help():
        help_text = """
        Usage: 
            python freeEnergyLandscape.py path/to/cv1_data.txt path/to/cv2_data.txt

        Optional arguments:
            --temperature           [int]       Simulation temperature in Kelvin (default: 300K)
            --kb                    [float]     Boltzmann constant in kJ/(mol·K) (default: 8.314e-3)
            --energy                [int]       Energy, single value (default: None)
            --bins_energy_histogram [int]       Bins for energy histogram (default: 100)
            --kde_bandwidth         [float]     Bandwidth for kernel density estimation (default: None)
            --names                 [str] [str] Names for the collective variables (default: CV1, CV2)
            --gif_angles            [int]       Angles for 3D GIF rotation (default: 10)
            --gif_elevation         [int]       Elevation angle for the 3D GIF (default: 10)
            --gif_duration          [float]     Duration per frame in the GIF in seconds (default: 0.1)

        Example:
            python freeEnergyLandscape.py cv1.txt cv2.txt --names Angle_CV1 Distance_CV2 --temperature 310 --energy 5 --bins_energy_histogram 100 --kde_bandwidth 0.5 --gif_angles 20

        """
        #      Notes:
        #  - The --energy argument can be a single value (e.g., 10) to indicate a maximum energy threshold, or a list of tuples to specify multiple energy ranges.
        #  - Use quotes for arguments that include spaces or special characters (e.g., --energy "[(0, 1), (1, 2)]").
        print(help_text)


    def main(self, energy_threshold, cv_names, n_angles, elevation, duration_per_frame):
        
        # Verificar ambos os arquivos de entrada antes de carregar os dados
        self.verify_input(self.cv1_path)
        self.verify_input(self.cv2_path)
        self.load_data()
        
        self.boltzmann_inversion(
            data_list=[self.proj1_data_original, self.proj2_data_original], 
            titles=cv_names, 
            threshold=energy_threshold
            )
        
        self.plot_histogram(
            data_list=[self.proj1_data_original, self.proj2_data_original], 
            titles=cv_names
            )

        self.cv_by_frame(
            data_list=[self.proj1_data_original, self.proj2_data_original], 
            titles=cv_names
            )

        self.plot_energy_landscape(
            threshold=energy_threshold, titles=cv_names
            )

        self.plot_3D_energy_landscape(
            threshold=energy_threshold, titles=cv_names
                                      )
        
        self.create_3D_gif(
            n_angles=n_angles, elevation=elevation, 
            duration_per_frame=duration_per_frame,
            titles=cv_names
                           )
        
        # self.save_low_energy_points_to_tsv(threshold=energy_threshold) # save low energy frames to tsv

        # Após o uso final dos dados, limpe-os para liberar memória
        if hasattr(self, 'cached_results'):
            del self.cached_results


if __name__ == "__main__":

    # Definindo valores padrão
    t = 300                     # --temperature           [int] [Kelvin]
    kB = 8.314e-3               # --kb                    [float] [kJ/(mol·K)]
    energy_threshold = None     # --energy                [int] [kJ/mol]
    bins_energy_histogram = 100 # --bins_energy_histogram [int]
    kde_bandwidth_cv = None     # --kde_bandwidth         [float]
    cv_names = ['CV1', 'CV2']   # --name                  [str] [str]
    n_angles = 10               # --gif_angles            [int]
    elevation = 10              # --gif_elevation         [int]
    duration_per_frame = 0.1    # --gif_duration          [float]

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
            elif key == "--bins_energy_histogram":
                bins_energy_histogram = int(sys.argv[i + 1])
                i += 2
            elif key == "--kde_bandwidth":
                kde_bandwidth_cv = float(sys.argv[i + 1]) if sys.argv[i + 1].lower() != "none" else None
                i += 2
            elif key == "--names":
                cv_names = [sys.argv[i + 1], sys.argv[i + 2]]  # Captura os nomes das variáveis coletivas
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
            else:
                print(f"Unrecognized option: {key}")
                sys.exit(1)
    else:
        FreeEnergyLandscape.help()
        sys.exit(1)

    try:
        fel = FreeEnergyLandscape(cv1_path, cv2_path, t, kB,  
                                bins=bins_energy_histogram, 
                                kde_bandwidth=kde_bandwidth_cv)
        
        fel.main(energy_threshold, cv_names=cv_names, 
                 n_angles=n_angles, elevation=elevation, 
                duration_per_frame=duration_per_frame)

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
