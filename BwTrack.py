
import pyBigWig
import sys

import numpy as np

class BwTrack:
    def __init__(self, bw_pl_path, bw_mn_path, signle_bw=False):
        '''
        Initialize the BwTrack object.
        
        Keyword arguments:
        - bw_pl_path: path to the plus strand bigwig file
        - bw_mn_path: path to the minus strand bigwig file
        - single_bw: if there is only plus strand bigwig file
        '''
        self.single_bw = signle_bw

        self.bw_pl = pyBigWig.open(bw_pl_path)

        if not self.single_bw:
            self.bw_mn = pyBigWig.open(bw_mn_path)
        else:
            self.bw_mn = self.bw_pl # for compatibility
        
    def __del__(self):
        '''
        Close the bigwig files.
        '''
        self.bw_pl.close()
        if not self.single_bw:
            self.bw_mn.close()

    @staticmethod
    def quantify_signal(signal: np.array, output_type: str):
        '''
        Quantify the signal. For quantification that is 
        not mirror-invariant, the signal should run from 
        the initiation site to the termination site.

        Keyword arguments:
        - signal: signal to quantify
        - output_type: what information is outputted (see BwTrack.get_supported_quantification_type)
        '''
        if output_type == "raw_count":
            return np.sum(signal)
        elif output_type == "RPK":
            return np.sum(signal) / len(signal) * 1e3
        else:
            raise Exception("Unsupported output type ({}).".format(output_type))
    
    @staticmethod
    def get_supported_quantification_type():
        '''
        Return the list of quantification types.
        '''
        return ["raw_count", "RPK"]


    def count_single_region(self, chrom, start, end, strand, 
                            output_type="raw_count", 
                            l_pad=0, r_pad=0, 
                            min_len_after_padding=50,
                            method_resolving_invalid_padding="raise",
                            ):
        '''
        Count the reads in a single region.
        Return the counts.

        Keyword arguments:
        - chrom: chromosome
        - start: start position
        - end: end position
        - strand: strandness, "+" or "-" or "."
        - output_type: what information is outputted [raw_count, RPK]
        - l_pad: left padding
        - r_pad: right padding
        - min_len_after_padding: minimum length of the region after padding
        - method_resolving_invalid_padding: method to resolve invalid padding
        '''
        # Invariants
        if self.single_bw:
            # For single_bw, no strand specificity would be possible
            assert strand == "."

        # Processing Padding
        if - l_pad - r_pad + min_len_after_padding > end - start:
            if method_resolving_invalid_padding == "fallback":
                start = start
                end = end
            elif method_resolving_invalid_padding == "raise":
                raise Exception("Padding is larger than the region ({}:{:d}-{:d}).".format(chrom, start, end))
            elif method_resolving_invalid_padding == "drop":
                sys.stderr.write("Padding is larger than the region ({}:{:d}-{:d}). Dropping the region.\n".format(chrom, start, end))
                return np.nan
            else:
                raise Exception("Unsupported method to resolve invalid padding ({}).".format(method_resolving_invalid_padding))
        else: 
            start = start - l_pad
            end = end + r_pad

        # Count signal
        if strand == "+":
            pl_sig = np.nan_to_num(self.bw_pl.values(chrom, 
                                                     start, 
                                                     end, 
                                                     ))
            quantification = BwTrack.quantify_signal(pl_sig, 
                                                     output_type, 
                                                     )
        elif strand == "-":
            mn_sig = -np.nan_to_num(self.bw_mn.values(chrom, 
                                                      start, 
                                                      end, 
                                                      ))
        
            quantification = BwTrack.quantify_signal(mn_sig, 
                                                     output_type, 
                                                     )
        elif strand == ".":
            pl_sig = np.nan_to_num(self.bw_pl.values(chrom, 
                                                     start, 
                                                     end))
            
            quantification = BwTrack.quantify_signal(pl_sig, 
                                                     output_type, 
                                                     )

            if not self.single_bw:
                mn_sig = -np.nan_to_num(self.bw_mn.values(chrom, 
                                                          start, 
                                                          end))
                mn_sig = np.flip(mn_sig)
                
                quantification += BwTrack.quantify_signal(mn_sig, 
                                                          output_type, 
                                                          )
        else:
            raise Exception("Invalid strand type.")

        return quantification

