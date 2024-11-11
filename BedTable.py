
# bed table classes to store bed files as pd.DataFrame

import sys

import numpy as np
import pandas as pd

from .exceptions import BedTableLoadException

# Only pd.Int64Dtype(), np.float64, and "O" are supported in this class
TYPE2PDDTYPE = {int: pd.Int64Dtype(), 
                np.int32: pd.Int64Dtype(),
                np.int64: pd.Int64Dtype(),
                pd.Int64Dtype(): pd.Int64Dtype(),
                float: np.float64, 
                np.float32: np.float64,
                np.float64: np.float64,
                str: pd.StringDtype(), 
                pd.StringDtype(): pd.StringDtype(),
                }

PDDTYPE2NANVAL = {pd.Int64Dtype(): np.nan,
                  np.float64: np.nan,
                  pd.StringDtype(): None,
                  }

class BedTableIterator:
    def __init__(self, bed_table: 'BedTable3') -> None:
        self.__bed_table = bed_table
        self.__index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.__index >= len(self.__bed_table):
            raise StopIteration
        else:
            region = self.__bed_table.get_region_by_index(self.__index)
            self.__index += 1
            return region

class BedRegion:
    '''
    A region type for bed file entires. 
    
    Enable easy access to bed file entries.
    - dict like access of chrom, start, end and other fields
    - support for comparison operations
    '''
    def __init__(self, chrom, start, end, **other_fields):
        self.__data_dict = {"chrom": chrom, "start": start, "end": end}
        self.__data_dict.update(other_fields)

    def to_dict(self):
        return self.__data_dict.copy()

    def __getitem__(self, key):
        return self.__data_dict[key]

    def __lt__(self, other):
        if self["chrom"] < other["chrom"]:
            return True
        elif self["chrom"] > other["chrom"]:
            return False

        # self["chrom"] == other["chrom"] case
        if self["start"] < other["start"]:
            return True
        elif self["start"] > other["start"]:
            return False
        
        # chrom and start are all the same
        return self["end"] < other["end"]
    
    def __eq__(self, other):
        return (self["chrom"] == other["chrom"]) and (self["start"] == other["start"]) and (self["end"] == other["end"])
    
    def __ne__(self, other):
        return not (self == other)
        
    def __gt__(self, other):
        return not ((self < other) or (self == other))
    
    def __le__(self, other):
        return (self < other) or (self == other)
    
    def __gt__(self, other):
        return not (self < other) and not (self == other)

class BedTable3:
    def __init__(self):
        self._data_df = pd.DataFrame(columns=self.column_names)

    # public methods
    @property
    def column_names(self):
        return ["chrom", "start", "end"]
    
    @property
    def column_types(self):
        '''
        A dictionary of column types.
        '''
        dtype_dict = {"chrom": str, "start": int, "end": int}
        
        dtype_dict = self._dtype2pddtype(dtype_dict)
        return dtype_dict
    
    @property
    def extra_column_names(self):
        return []
    
    @property
    def extra_column_dtype(self):
        return []

    def load_from_file(self, ipath: str) -> None:
        '''
        Load a bed file.
        '''
        try:
            self._data_df = pd.read_csv(ipath, 
                                        sep="\t", 
                                        names=self.column_names,
                                        na_values='.', 
                                        )
        except ValueError as e:
            raise BedTableLoadException(f"Error loading bed file: number of columns does not match.")
        
        self._force_dtype()

        if not self._is_sorted():
            self._sort()

    def load_from_dataframe(self, df: pd.DataFrame, 
                            column_map=None) -> None:
        '''
        Load a pd.DataFrame.

        Keyword arguments:
        - df: a pd.DataFrame that contains the data for a bed table
        - column_map: a dictionary that maps the column names in the 
                      bed table to the column names in the pd.DataFrame         
        '''
        if not column_map:
            column_map = {col: col for col in self.column_names}

        try:
            self._data_df = pd.DataFrame(df[[column_map[col] for col in self.column_names]].values, 
                                         columns=self.column_names, 
                                         )

        except ValueError as e:
            raise BedTableLoadException(f"Error loading pd.DataFrame: number of columns does not match.")
        
        self._force_dtype()

        if not self._is_sorted():
            self._sort()

    def apply_logical_filter(self, logical_array: np.array) -> 'BedTable3':
        '''
        Use a logical np.array to filter the table.
        Return a new BedTable3 instance.
        '''
        # input check removed for compatibility with newer pandas version
        #TODO: add updated input check
        # if not isinstance(logical_array, np.ndarray):
        #     raise ValueError("logical_array must be a np.array")
        # if not logical_array.dtype == bool:
        #     raise ValueError("logical_array must be a boolean np.array")

        new_bed_table = self._clone_empty()
        new_bed_table.load_from_dataframe(self._data_df.loc[logical_array])

        return new_bed_table
        
    def _get_subset_data_df(self, chrom: str, start: int, end: int) -> pd.DataFrame:
        '''
        Subset the data df to the given region.
        Return a new pd.DataFrame.
        '''
        #TODO: Add only subsetting by chrome option
        subset_data_df = self._data_df.loc[self._data_df["chrom"] == chrom]
        subset_data_df = subset_data_df.loc[subset_data_df["start"] >= start]
        subset_data_df = subset_data_df.loc[subset_data_df["end"] <= end]

        return subset_data_df.copy()

    def region_subset(self, chrom: str, start: int, end: int) -> 'BedTable3':
        '''
        Subset the table to the given region.
        Only return regions that are fully contained in the given region.
        Return a new BedTable3 instance.
        '''
        subset_data_df = self._get_subset_data_df(chrom, start, end)

        new_bed_table = self._clone_empty()
        new_bed_table.load_from_dataframe(subset_data_df)

        return new_bed_table

    def to_dataframe(self) -> pd.DataFrame:
        '''
        Return the table as a copy of pd.DataFrame
        '''
        return self._data_df.copy()

    def write(self, opath: str) -> None:
        '''
        Write the table to a bed file.
        '''
        df2write = self._data_df.copy()

        for col, col_dtype in self.column_types.items():
            if col_dtype == pd.StringDtype():
                df2write[col] = df2write[col].fillna(".")
            else:
                df2write[col] = df2write[col].astype(pd.StringDtype())
                df2write[col] = df2write[col].fillna(".")

        if opath == "stdout":
            opath = sys.stdout

        df2write.to_csv(opath, 
                        sep="\t", 
                        header=False, 
                        index=False, 
                        )

    def get_chrom_names(self) -> np.array:
        '''
        Return a np.array of chrom names.
        '''
        return self._data_df["chrom"].values

    def get_start_locs(self) -> np.array:
        '''
        Return a np.array of start locations.
        '''
        return self._data_df["start"].values

    def get_end_locs(self) -> np.array:
        '''
        Return a np.array of end locations.
        '''
        return self._data_df["end"].values
    
    def get_region_by_index(self, index: int) -> pd.Series:
        '''
        Return a region by index.
        '''
        region_series = self._data_df.iloc[index]

        return BedRegion(**region_series.to_dict()
                         )
    
    def iter_regions(self) -> tuple:
        '''
        Return a iterator of regions.
        '''
        return BedTableIterator(self)
    
    def search_region(self, chrom: str, start: int, end: int, overlapping_base=1) -> np.array:
        '''
        Search for a region and return the region index as an array.
        Return empty array if no data is found for the given region.

        Keyword arguments:
        - chrom: chromosome
        - start: start location
        - end: end location
        - overlapping_base: the number of overlapping bases required for a match
        '''
        # first narrow down with chrom
        chrom_left_bound, chrom_right_bound = self._search_array(self.get_chrom_names(), chrom)
        chrom_left_bound += 1

        search_bt = self.subset_by_index(np.arange(chrom_left_bound, chrom_right_bound))

        # then narrow down with end
        end_left_bound, end_right_bound = self._search_array(search_bt.get_start_locs(), end)

        search_space = [chrom_left_bound, chrom_left_bound + end_left_bound + 1]

        result_inds = []
        for i in range(search_space[0], search_space[1]):
            region = self.get_region_by_index(i)

            if region["end"] - start >= overlapping_base and end - region["start"] >= overlapping_base:
                result_inds.append(i)
        
        return np.array(result_inds, dtype=int)
    
    def concat(self, other_bed_table):
        '''
        Concatenate two bed tables.
        Return a new BedTable instance.
        '''
        new_bed_table = self._clone_empty()
        new_bed_table.load_from_dataframe(pd.concat([self._data_df, other_bed_table._data_df], ignore_index=True))
        #TODO: Use merge sort to improve performance

        return new_bed_table
    
    def subset_by_index(self, index_array: np.array):
        '''
        Subset the table by index array.
        Return a new BedTable instance.
        '''
        new_bed_table = self._clone_empty()
        new_bed_table.load_from_dataframe(self._data_df.loc[index_array].copy())

        return new_bed_table
    
    def _clone_empty(self):
        '''
        Clone an empty instance.
        '''
        new_bed_table = self.__class__()
        return new_bed_table

    def _sort(self) -> None:
        '''
        Sort the bed table.
        '''
        self._data_df.sort_values(by=["chrom", "start", "end"], inplace=True)
        self._data_df.reset_index(drop=True, inplace=True)

    def _is_sorted(self) -> bool:
        '''
        Check if the table is sorted.
        '''
        for i in range(1, len(self)):
            if self.get_region_by_index(i) < self.get_region_by_index(i-1):
                return False
        
        return True

    def _dtype2pddtype(self, dtype_dict):
        '''
        Convert dtype to pd.Dtype.
        '''
        for key in dtype_dict.keys():
            if dtype_dict[key] in TYPE2PDDTYPE.keys():
                dtype_dict[key] = TYPE2PDDTYPE[dtype_dict[key]]
            else:
                raise ValueError(f"Unsupported dtype: {dtype_dict[key]}")
        return dtype_dict

    def _force_dtype(self):
        '''
        Force the column types.
        '''
        for field, field_dtype in self.column_types.items():
            self._data_df[field] = self._data_df[field].replace(".", PDDTYPE2NANVAL[field_dtype])
            self._data_df[field] = self._data_df[field].astype(field_dtype)

    def __len__(self) -> int:
        return self._data_df.shape[0]
    
    @staticmethod
    def _search_array(array: np.array, value) -> np.array:
        '''
        Search for a value in an sorted array and return 2 indexes.
        ind1 will be the index for last element smaller than value, 
        and ind2 will be the index for first element larger than value.

        If the first element is larger than value, ind1 will be -1.
        If the last element is smaller than value, ind2 will be len(array).

        Keyword arguments:
        - array: a np.array
        - value: the value to search

        Return a tuple of two integers.
        '''
        # binary search for ind 1
        search_start = -1
        search_end = len(array)

        while search_start+1 < search_end:
            new_ind = (search_start + search_end) // 2
            if array[new_ind] < value:
                search_start = new_ind
            else:
                search_end = new_ind
        
        ind1 = search_start

        # binary search for ind 2
        search_start = -1
        search_end = len(array)

        while search_start+1 < search_end:
            new_ind = (search_start + search_end) // 2
            if array[new_ind] <= value:
                search_start = new_ind
            else:
                search_end = new_ind
        
        ind2 = search_end

        return ind1, ind2


class BedTable6(BedTable3):
    def __init__(self):
        self._data_df = pd.DataFrame(columns=self.column_names)

    @property
    def column_names(self):
        return ["chrom", "start", "end", "name", "score", "strand"]

    @property
    def column_types(self):
        column_type = super().column_types
        column_type["name"] = str
        column_type["score"] = float
        column_type["strand"] = str

        self._dtype2pddtype(column_type)
        return column_type
    
    def get_region_names(self) -> np.array:
        '''
        Return a np.array of region names.
        '''
        return self._data_df["name"].values
    
    def get_region_scores(self) -> np.array:
        '''
        Return a np.array of region scores.
        '''
        return self._data_df["score"].values
    
    def get_region_strands(self) -> np.array:
        '''
        Return a np.array of region strands.
        '''
        return self._data_df["strand"].values
    
    def region_subset(self, chrom: str, start: int, end: int) -> 'BedTable6':
        subset_data_df = self._get_subset_data_df(chrom, start, end)

        new_bed_table = self._clone_empty()
        new_bed_table.load_from_dataframe(subset_data_df)

        return new_bed_table
    
    def load_from_BedTable3(self, bed3: BedTable3) -> None:
        '''
        Convert a BedTable3 to BedTable6.
        '''
        self._data_df = bed3.to_dataframe()
        self._data_df["name"] = "."
        self._data_df["score"] = "."
        self._data_df["strand"] = "."

        self._force_dtype()
    
class BedTable6Plus(BedTable6):
    def __init__(self, 
                 extra_column_names: list,
                 extra_column_dtype: list = None,
                 ):
        self._extra_column_names = extra_column_names

        if not extra_column_dtype:
            self._extra_column_dtype = [str] * len(extra_column_names)
        else:
            self._extra_column_dtype = extra_column_dtype
        
        self._data_df = pd.DataFrame(columns=self.column_names)

    @property
    def column_names(self):
        return ["chrom", "start", "end", "name", "score", "strand"] + self.extra_column_names
    
    @property
    def column_types(self):
        column_type = super().column_types
        for extra_col, extra_col_dtype in zip(self.extra_column_names, self.extra_column_dtype):
            column_type[extra_col] = extra_col_dtype
        
        self._dtype2pddtype(column_type)

        return column_type
    
    @property
    def extra_column_names(self):
        return self._extra_column_names
    
    @property
    def extra_column_dtype(self):
        return self._extra_column_dtype

    def get_region_extra_column(self, column_name) -> np.array:
        '''
        Return a np.array of extra column data for all the regions. Given the column name.
        '''
        return self._data_df[column_name].values
    
    def _clone_empty(self):
        new_bed_table = self.__class__(self.extra_column_names, self.extra_column_dtype)
        return new_bed_table

class BedTablePairEnd(BedTable3):
    def __init__(self, 
                 extra_column_names: list = None,
                 extra_column_dtype: list = None,
                 ):
        if not extra_column_names:
            self._extra_column_names = []
        else:
            self._extra_column_names = extra_column_names

        if not extra_column_dtype:
            self._extra_column_dtype = [str] * len(extra_column_names)
        else:
            self._extra_column_dtype = extra_column_dtype
        
        self._data_df = pd.DataFrame(columns=self.column_names)
        self._other_region_inverse_index = self._build_inverse_index_for_the_other_region()
    
    def load_from_dataframe(self, df: pd.DataFrame, column_map=None) -> None:
        super().load_from_dataframe(df, column_map)
        self._other_region_inverse_index = self._build_inverse_index_for_the_other_region()

    def load_from_file(self, ipath: str) -> None:
        super().load_from_file(ipath)
        self._other_region_inverse_index = self._build_inverse_index_for_the_other_region()

    @property
    def extra_column_names(self):
        return self._extra_column_names
    
    @property
    def extra_column_dtype(self):
        return self._extra_column_dtype

    @property
    def column_names(self):
        return ["chrom", "start", "end", "chrom2", "start2", "end2", "name", "score", "strand", "strand2"] + self.extra_column_names
    
    @property
    def column_types(self):
        column_type = super().column_types
        column_type["chrom2"] = str
        column_type["start2"] = int
        column_type["end2"] = int
        column_type["name"] = str
        column_type["score"] = float
        column_type["strand"] = str
        column_type["strand2"] = str

        for extra_col, extra_col_dtype in zip(self.extra_column_names, self.extra_column_dtype):
            column_type[extra_col] = extra_col_dtype

        self._dtype2pddtype(column_type)

        return column_type
    
    def get_other_region_chroms(self) -> np.array:
        '''
        Return a np.array of chroms of the other region.
        '''
        return self._data_df["chrom2"].values
    
    def get_other_region_starts(self) -> np.array:
        '''
        Return a np.array of starts of the other region.
        '''
        return self._data_df["start2"].values
    
    def get_other_region_ends(self) -> np.array:
        '''
        Return a np.array of ends of the other region.
        '''
        return self._data_df["end2"].values
    
    def get_pair_names(self) -> np.array:
        '''
        Return a np.array of region names.
        '''
        return self._data_df["name"].values
    
    def get_pair_scores(self) -> np.array:
        '''
        Return a np.array of pair end region scores.
        '''
        return self._data_df["score"].values
    
    def get_region_strands(self) -> np.array:
        '''
        Return a np.array of region strands.
        '''
        return self._data_df["strand"].values
    
    def get_other_region_strands(self) -> np.array:
        '''
        Return a np.array of strands of the other region in the pair.
        '''
        return self._data_df["strand2"].values

    def get_extra_column(self, column_name) -> np.array:
        '''
        Return a np.array of extra column data for all the regions. Given the column name.
        '''
        return self._data_df[column_name].values
    
    def search_pair_extra_column(self, chr1, start1, end1, 
                                 chr2, start2, end2, 
                                 column_name, 
                                 overlapping_base1=1, 
                                 overlapping_base2=1, 
                                 ):
        '''
        Search for a pair end region and return the column data as indicated by the column name.
        Return None if no data is found for the given pair regions.

        Keyword arguments:
        - chr1: chromosome of the first region
        - start1: start of the first region
        - end1: end of the first region
        - chr2: chromosome of the second region
        - start2: start of the second region
        - end2: end of the second region
        - column_name: the column name of the data to return
        - overlapping_base1: the number of overlapping bases required for a match for the first input region
        - overlapping_base2: the number of overlapping bases required for a match for the second input region
        '''
        # forward search
        temp_ind_array = self.search_region(chr1, start1, end1, overlapping_base1)
        temp_bed_table = self.subset_by_index(temp_ind_array)
        forward_inds = temp_bed_table.search_second_region(chr2, start2, end2, overlapping_base2)
        forward_result_bed_table = temp_bed_table.subset_by_index(forward_inds)

        # backward search
        temp_ind_array = self.search_region(chr2, start2, end2, overlapping_base2)
        temp_bed_table = self.subset_by_index(temp_ind_array)
        backward_inds = temp_bed_table.search_second_region(chr1, start1, end1, overlapping_base1)
        backward_result_bed_table = temp_bed_table.subset_by_index(backward_inds)

        all_result_bt = forward_result_bed_table.concat(backward_result_bed_table)

        return all_result_bt.get_extra_column(column_name)

    def search_second_region(self, chrom: str, start: int, end: int, overlapping_base=1) -> np.array:
        '''
        Search for the second region in the pair and return the index as an array.
        Return empty array if no data is found for the given region.

        Keyword arguments:
        - chrom: chromosome
        - start: start location
        - end: end location
        - overlapping_base: the number of overlapping bases required for a match
        '''
        searched_ind = self._other_region_inverse_index.search_region(chrom, start, end, overlapping_base)
        return self._other_region_inverse_index.get_region_extra_column("index")[searched_ind]

    def _clone_empty(self):
        new_bed_table = self.__class__(self.extra_column_names, self.extra_column_dtype)
        return new_bed_table

    def _build_inverse_index_for_the_other_region(self):
        '''
        Build an inverse index for the other region.
        '''
        inverse_index_df = self._data_df.loc[:, ["chrom2", "start2", "end2", "name", "score", "strand2"]].copy()
        inverse_index_df["index"] = self._data_df.index

        inverse_index_bt = BedTable6Plus(extra_column_names=["index"], 
                                         extra_column_dtype=[int], 
                                         )
        
        inverse_index_bt.load_from_dataframe(inverse_index_df, 
                                             column_map={"chrom": "chrom2", 
                                                         "start": "start2", 
                                                         "end": "end2", 
                                                         "name": "name", 
                                                         "score": "score", 
                                                         "strand": "strand2", 
                                                         "index": "index", 
                                                         }
                                             )

        return inverse_index_bt