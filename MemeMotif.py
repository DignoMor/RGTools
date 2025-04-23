
class MemeMotif:
    def __init__(self, file_path: str):
        '''
        Initializes the MemeMotif object with the given file path.

        Keyword arguments:
        - file_path: Path to the MEME motif file.
        '''
        pass

    def _parse_meme_file(self, file_path: str):
        '''
        Parses the MEME motif file and extracts relevant information.
        Part of the initialization process.
        
        Keyword arguments:
        - file_path: Path to the MEME motif file.
        '''
        pass

    def get_meme_version(self):
        '''
        Returns the version of the MEME motif file.

        Return:
        - A string representing the version of the MEME motif file.
        '''
        pass

    def get_alphabet(self):
        '''
        Returns the alphabet used in the MEME motif file.

        Return: 
        - A string representing the alphabet used in the MEME motif file.
        '''
        pass

    def get_strands(self):
        '''
        Returns the strands used in the MEME motif file.

        Return:
        - A list of strings representing the strands used in the MEME motif file.
        '''
        pass

    def get_bg_freq(self):
        '''
        Returns the background frequency used in the MEME motif file.

        Return: 
        - A dictionary with keys being the alphabet characters and values being their frequencies.
        '''
        pass

    def get_motif_list(self):
        '''
        Returns the list of motifs in the MEME motif file.

        Return: 
        - A list of strings representing the names of the motifs in the MEME motif file.
        '''
        pass

    def get_motif_pwm(self, motif_name):
        '''
        Returns the position weight matrix (PWM) of the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return:
        - A 2D numpy array representing the PWM of the given motif.
          The array should have shape (num_positions, num_alphabet_chars).
        '''
        pass

    def get_motif_alphabet_length(self, motif_name):
        '''
        Returns the length of the alphabet used in the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return: 
        - An integer representing the length of the alphabet used in the given motif.
        '''
        pass

    def get_motif_length(self, motif_name):
        '''
        Returns the length of the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return:
        - An integer representing the length of the given motif.
        '''
        pass

    def get_motif_num_source_sites(self, motif_name):
        '''
        Returns the number of source sites for the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return: 
        - An integer representing the number of source sites for the given motif.
        '''
        pass

    def get_motif_source_eval(self, motif_name):
        '''
        Returns the source evaluation for the given motif.

        Keyword arguments:
        - motif_name: Name of the motif.

        Return: 
        - A float representing the source E-value for the given motif.
        '''
        pass
    

