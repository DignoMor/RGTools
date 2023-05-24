# only works for GENCODE now

import pandas as pd

class GTFRecord:
    def __init__(self, chr_name:str, record_source:str, feature_type:str, start_loc:int, 
                 end_loc:int, score:str, strand:str, phase:str, add_info:dict):
        self.__general_info_dict = dict()
        
        self.__general_info_dict["chr_name"] = chr_name
        self.__general_info_dict["record_source"] = record_source
        self.__general_info_dict["feature_type"] = feature_type
        self.__general_info_dict["start_loc"] = int(start_loc)
        self.__general_info_dict["end_loc"] = int(end_loc)
        self.__general_info_dict["score"] = score
        self.__general_info_dict["strand"] = strand
        self.__general_info_dict["phase"] = phase

        self.__add_info_dict = add_info
    
        # sanity check
        # TODO: add more
        assert self.__general_info_dict["strand"] in ("+", "-")
    
    def _print_all_info(self):
        print(self.__general_info_dict)
        print(self.__add_info_dict)

    def search_general_info(self, query):
        if query in self.__general_info_dict.keys():
            return self.__general_info_dict[query]
        else:
            return None

    def search_add_info(self, query):
        if query in self.__add_info_dict.keys():
            return self.__add_info_dict[query]
        else:
            return None
    
class GTFHandle:
    def __init__(self, gtf_path, filter=lambda x: True):
        self.__record_filter = filter
        self.gtf_path = gtf_path
        self.gtf_file = open(gtf_path, 'r')

        self.__first_line = self.gtf_file.readline()

        self.comments=[]
        while self.__first_line[0] == "#":
            self.comments.append(self.__first_line)
            self.__first_line = self.gtf_file.readline()
        
        self.gtf_file.close()

    def __iter__(self):
        self.gtf_file = open(self.gtf_path, 'r')
        self.__first_line = self.gtf_file.readline()

        while self.__first_line[0] == "#":
            self.__first_line = self.gtf_file.readline()

        self.__current_line = self.__first_line
        return self

    def __next__(self):
        while True: 
            # skip empty lines
            while self.__current_line == "\n":
                self.__current_line = self.gtf_file.readline()
            
            # Stop iteration if eof
            if not self.__current_line: 
                self.gtf_file.close()
                raise StopIteration

            current_record = self.__parse_line(self.__current_line)
            self.__current_line = self.gtf_file.readline()
            
            if self._is_record_passed_filter(current_record):
                # start again if not passing the filter
                break
        
        return current_record

    def __parse_line(self, line_to_parse):
        chr_name, record_source, feature_type, start_loc, end_loc, \
        score, strand, phase, add_info = \
            line_to_parse.split("\t")
        
        add_info_dict = dict()

        for add_info_field in add_info.replace("\n", "").rstrip(";").split(";"):
            field_name = add_info_field.split()[0].strip('"')
            field_value = add_info_field.split()[1].strip('"')
            add_info_dict[field_name] = field_value
        
        current_record = GTFRecord(chr_name=chr_name, 
                                   record_source=record_source, 
                                   feature_type=feature_type, 
                                   start_loc=int(start_loc), 
                                   end_loc=int(end_loc), 
                                   score=score, 
                                   strand=strand, 
                                   phase=phase, 
                                   add_info=add_info_dict,
                                   )

        return current_record
    
    def _is_record_passed_filter(self, record):
        # Return False if to filter
        return self.__record_filter(record)
        
    def get_comments(self):
        comment_lines = [s[1:].lstrip().replace("\n", "") for s in self.comments]
        return "\n".join(comment_lines)
    
    def filter_by_general_record(self, key, value):
        new_record_filter = lambda r: r.search_general_info(key) == value
        record_filter = lambda r: (new_record_filter(r) & self.__record_filter(r))
        return GTFHandle(self.gtf_path, filter=record_filter)

    def filter_by_add_record(self, key, value):
        new_record_filter = lambda r: r.search_add_info(key) == value
        record_filter = lambda r: (new_record_filter(r) & self.__record_filter(r))
        return GTFHandle(self.gtf_path, filter=record_filter)
    
    def count_total(self):
        count = 0
        for record in self:
            count += 1
        
        return count
    
    def write_tss_bed(self, opath, l_padding=0, r_padding=0):
        tss_list = []
        gene_id_list = []
        chr_name_list = []
        score_list = []
        strand_list = []

        for record in self:
            if record.search_general_info("strand") == "+":
                tss = record.search_general_info("start_loc") 
            else: 
                tss = record.search_general_info("end_loc") - 1 # GTF is closed at start and open at end

            gene_id = record.search_add_info("gene_id")
            chr_name = record.search_general_info("chr_name")
            score = record.search_general_info("score")
            strand = record.search_general_info("strand")

            tss_list.append(tss)
            gene_id_list.append(gene_id)
            chr_name_list.append(chr_name)
            score_list.append(score)
            strand_list.append(strand)
                
        bed_df = pd.DataFrame({
            "chr_name": chr_name_list, 
            "start": [t - l_padding - 1 for t in tss_list], # bed is 0-based
            "end": [t + r_padding for t in tss_list], 
            "gene_id": gene_id_list, 
            "score": score_list, 
            "strand": strand_list, 
        })

        bed_df.to_csv(opath, 
                      sep="\t", 
                      header=False, 
                      index=False, 
                      )

