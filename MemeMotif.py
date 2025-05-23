
import re 
import numpy as np

class MemeMotif:
    def __init__(self, file_path: str):
        '''
        Initializes the MemeMotif object with the given file path.

        Keyword arguments:
        - file_path: Path to the MEME motif file.
        '''
        self.file_path = file_path
        self.version = None
        self.alphabet = None
        self.strands = None
        self.bg_freq = None
        self.motifs = []
        self.motif_info_dict = {}

        self._parse_meme_file(file_path)
        

    def _parse_meme_file(self, file_path: str):
        '''
        Parses the MEME motif file and extracts relevant information.
        Part of the initialization process.
        
        Keyword arguments:
        - file_path: Path to the MEME motif file.
        '''

        with open(file_path, "r") as f:
            while line := f.readline():
                line = line.strip()
                if line.startswith("MEME version"):
                    self.version = line.split()[-1].strip()
                    break
            
            while line := f.readline():
                line = line.strip()
                if line.startswith("ALPHABET"):
                    self.alphabet = line.split("=")[1].strip()
                    break

            while line := f.readline():
                line = line.strip()
                if line.startswith("strands"):
                    self.strands = line.split()[1:]
                    break

            while line := f.readline():
                line = line.strip()
                if line.startswith("Background letter frequencies"):
                    line = f.readline().strip()
                    line_info = line.split()

                    freq_dict = {}

                    for i in range(len(self.alphabet)):
                        char = line_info[2*i]
                        freq = line_info[2*i + 1]

                        freq_dict[char] = freq
                    
                    self.bg_freq = [float(freq_dict[c]) for c in self.alphabet]
                    break

            while line := f.readline():
                line = line.strip()
                if line.startswith("MOTIF"):
                    motif_name = line.split()[1]
                    motif_info = {}

                    line = f.readline().strip()
                    match = re.match(r"letter-probability matrix: alength=(.*)w=(.*)nsites=(.*)E=(.*)", line)

                    motif_info["alphabet_length"] = int(match.group(1))
                    motif_info["motif_length"] = int(match.group(2))
                    motif_info["num_source_sites"] = int(match.group(3))
                    motif_info["source_eval"] = float(match.group(4))

                    pwm_arr = np.zeros((motif_info["motif_length"], motif_info["alphabet_length"]))
                    for i in range(motif_info["motif_length"]):
                        line = f.readline().strip()
                        line_info = line.split()
                        pwm_arr[i] = np.array(line_info, dtype=float)
                
                    motif_info["pwm"] = pwm_arr

                    self.motifs.append(motif_name)
                    self.motif_info_dict[motif_name] = motif_info

    def get_meme_version(self):
        '''
        Returns the version of the MEME motif file.

        Return:
        - A string representing the version of the MEME motif file.
        '''
        return self.version

    def get_alphabet(self):
        '''
        Returns the alphabet used in the MEME motif file.

        Return: 
        - A string representing the alphabet used in the MEME motif file.
        '''
        return self.alphabet

    def get_strands(self):
        '''
        Returns the strands used in the MEME motif file.

        Return:
        - A list of strings representing the strands used in the MEME motif file.
        '''
        return self.strands

    def get_bg_freq(self):
        '''
        Returns the background frequency used in the MEME motif file.

        Return: 
        - A dictionary with keys being the alphabet characters and values being their frequencies.
        '''
        return self.bg_freq

    def get_motif_list(self):
        '''
        Returns the list of motifs in the MEME motif file.

        Return: 
        - A list of strings representing the names of the motifs in the MEME motif file.
        '''
        return self.motifs

    def get_motif_pwm(self, motif_name):
        '''
        Returns the position weight matrix (PWM) of the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return:
        - A 2D numpy array representing the PWM of the given motif.
          The array should have shape (num_positions, num_alphabet_chars).
        '''
        return self.motif_info_dict[motif_name]["pwm"]

    def get_motif_alphabet_length(self, motif_name):
        '''
        Returns the length of the alphabet used in the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return: 
        - An integer representing the length of the alphabet used in the given motif.
        '''
        return self.motif_info_dict[motif_name]["alphabet_length"]

    def get_motif_length(self, motif_name):
        '''
        Returns the length of the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return:
        - An integer representing the length of the given motif.
        '''
        return self.motif_info_dict[motif_name]["motif_length"]

    def get_motif_num_source_sites(self, motif_name):
        '''
        Returns the number of source sites for the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return: 
        - An integer representing the number of source sites for the given motif.
        '''
        return self.motif_info_dict[motif_name]["num_source_sites"] 

    def get_motif_source_eval(self, motif_name):
        '''
        Returns the source evaluation for the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return: 
        - A float representing the source E-value for the given motif.
        '''
        return self.motif_info_dict[motif_name]["source_eval"]
    

