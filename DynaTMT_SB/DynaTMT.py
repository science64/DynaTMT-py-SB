'''

    DynaTMT-py - a python package to process SILAC/TMT proteomics data
    Copyright (C) 2021  Kevin Klann - 2023 Süleyman Bozkurt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''

__author__ = "Kevin Klann - Süleyman Bozkurt"
__version__ = "v2.7.0"
__maintainer__ = "Süleyman Bozkurt"
__email__ = "sbozkurt.mbg@gmail.com"
__date__ = '18.01.2021'
__update__ = '05.12.2023'

from scipy.stats import trim_mean
import pandas as pd
import numpy as np
from numpy.random import random
import warnings
warnings.filterwarnings("ignore")

class PD_input:
    '''Class containing functions to analyze the default peptide/PSM output from ProteomeDiscoverer. All column names are assumed
    to be default PD output. If your column names do not match these assumed strings, you can modify them or use the plain_text_input
    class, that uses column order instead of names.
    '''
    def __init__(self, input):
        '''Initialises PD_input class with specified input file. The input file gets stored
        in the class variable self.input_file.
        '''
        # self.input_file = input # this has been changed, we are not using it anymore!
        self.channels = self.get_channels(input)

    def get_channels(self, input):
        channels = [col for col in input.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in input.columns if 'Abundance' in col]
        return channels

    def log_func(func):  # to show (log) what function is under usage!
        def wrapper(self, *args, **kwargs):
            print(f"Calling function: {func.__name__}")
            return func(self, *args, **kwargs)
        return wrapper

    @log_func
    def filter_peptides(self, filtered_input):
        filtered_input = filtered_input[~filtered_input['Master Protein Accessions'].str.contains(';',na=False)]
        filtered_input = filtered_input[filtered_input['Contaminant'] == False]

        # this part removes the peptides with 0 signal/noise ratio (meaning that all the channels are 0)
        try:
            filtered_input = filtered_input.dropna(subset=['Average Reporter SN'])
            filtered_input = filtered_input[filtered_input['Average Reporter SN'] != 0]
        except:
            filtered_input = filtered_input.dropna(subset=['Average Reporter S/N'])
            filtered_input = filtered_input[filtered_input['Average Reporter S/N'] != 0]

        # this part removoes whole channel sum of 0, because we remove booster and rest might be 0 but Average Reporter SN is not 0
        # Calculate the row-wise sum for the specified columns
        row_sums = filtered_input[self.channels].sum(axis=1)
        # Filter out rows where the sum is zero
        filtered_input = filtered_input[row_sums != 0]

        # this part replace NaN values with 0
        filtered_input.fillna(0, inplace=True)

        return filtered_input

    @log_func
    def IT_adjustment(self, input):
        '''This function adjusts the input DataFrame stored in the class variable self.input_file for Ion injection times.
        Abundance channels should contain "Abundance:" string and injection time uses "Ion Inject Time" as used by ProteomeDiscoverer
        default output. For other column headers please refer to plain_text_input class.
        '''
        IT=[col for col in input.columns if 'Ion Inject Time' in col]
        inject_times=input[IT[0]]
        input[self.channels]=input[self.channels].divide(inject_times,axis=0)
        input[self.channels]=input[self.channels].multiply(1000)

        print("IT adjustment done!")
        return input

    @log_func
    def total_intensity_normalisation(self, input):
        '''This function normalizes the self.input_file variable to the summed intensity of all TMT channels. It modifies the self.input_file
        to the updated DataFrame containing the normalized values.
        '''
        input_df = input.copy(deep=True)
        input_df = input_df.dropna(subset=self.channels)
        minimum = np.argmin(input_df[self.channels].sum().values)
        summed = np.array(input_df[self.channels].sum().values)
        minimum = summed[minimum]
        norm_factors = summed / minimum
        input_df.loc[:, self.channels] = input_df[self.channels].divide(norm_factors, axis=1)
        print("Total intensity normalisation done!")
        return input_df

    @log_func
    def Median_normalisation(self, input):
        '''This function normalizes the self.input_file variable to the median of all individual TMT channels. It modifies the self.input_file
        to the updated DataFrame containing the normalized values.
        '''
        input=input.dropna(subset=self.channels)
        minimum=np.argmin(input[self.channels].median().values)
        summed=np.array(input[self.channels].median().values)
        minimum=summed[minimum]
        norm_factors=summed/minimum
        input[self.channels]=input[self.channels].divide(norm_factors, axis=1)
        print("Median normalisation done!")
        return input

    @log_func
    def TMM(self, input):
        '''This function implements TMM normalisation (Robinson & Oshlack, 2010, Genome Biology). It modifies the self.input_file class
        variable.
        '''
        input=input.dropna(subset=self.channels)
        input_trim=input[input[self.channels] < input[self.channels].quantile(.95)]
        print("Normalization")
        input_trim[self.channels]=input_trim[self.channels].divide(input_trim[self.channels[0]],axis=0)
        tm=np.argmin(trim_mean(input_trim[self.channels],0.25))
        summed=np.array(trim_mean(input_trim[self.channels], 0.25))
        minimum=summed[tm]
        norm_factors=summed/minimum
        input[self.channels]=input[self.channels].divide(norm_factors, axis=1)
        return input

    @log_func
    def extract_heavy (self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all heavy labelled peptides. Naming of the 
        Modifications: Arg10: should contain Label, TMTK8, TMTproK8. Strings for modifications can be edited below for customisation.
        Returns heavy peptide DF
        '''
        modi = list([col for col in input.columns if 'Modification' in col])[0]
        '''Change Modification String here'''
        Heavy_peptides = input[input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6|TMTproK4|TMTproK6', na=False)]
        print("Extraction Done","Extracted Heavy Peptides:", len(Heavy_peptides))

        return Heavy_peptides

    @log_func
    def extract_light (self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all light labelled peptides. Naming of the 
        Modifications: Arg10: should contain Label, TMTK8, TMTproK8. Strings for modifications can be edited below for customisation.
        Returns light peptide DF
        '''
        modi = list([col for col in input.columns if 'Modification' in col])[0]
        light_peptides = input[~input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6',na=False)]
        print("Extraction Done","Extracted Light Peptides:", len(light_peptides))
        return light_peptides

    # @log_func
    # def baseline_correction(self, input_file, random=True, threshold=5, i_baseline=0, include_negatives=False):
    #     '''This function takes the input_file DataFrame and substracts the baseline/noise channel from all other samples. The index of the
    #     baseline column is defaulted to 0. Set i_baseline=X to change baseline column.
    #
    #     Threshold: After baseline substraction the remaining average signal has to be above threshold to be included. Parameter is set with threshold=X.
    #     This prevents very low remaining signal peptides to produce artificially high fold changes. Has to be determined empirically.
    #
    #     First PSMs combined into peptides by 'Master Protein Accessions', 'Annotated Sequence', 'Modifications' and then baseline correction done, export
    #     peptide file
    #     '''
    #
    #     random_float = np.random.RandomState(69)  # random seed for NaN, empty or 0 values.
    #
    #     # channels and 'Modifications', 'Master Protein Accessions', 'Annotated Sequence' are required for further analysis for baseline correction.
    #     peptide = input_file[self.channels + ['Modifications', 'Master Protein Accessions', 'Annotated Sequence']]
    #
    #     baseline_channel = self.channels[i_baseline]
    #     baseline = peptide[baseline_channel]
    #     peptide[self.channels] = peptide[self.channels].subtract(baseline, axis='index')
    #     peptide['Mean'] = peptide[self.channels].mean(axis=1)
    #
    #     # !!!!! this part might change !!!!
    #     peptide = peptide.loc[peptide['Mean'] >= threshold]  # set S/N threshold for each PSM
    #     peptide = peptide.drop("Mean", axis=1)
    #
    #     if (include_negatives == False and random == False):
    #         peptide[peptide < 0] = 0  # replace negative abundances with 0
    #
    #     elif (include_negatives == False and random == True):
    #         for channel in self.channels:
    #             peptide[channel] = np.where(peptide[channel] < 0, random_float.random_sample(size=len(peptide)),
    #                                         peptide[channel])
    #     else:  # for other conditions we are not doing anything.
    #         pass
    #
    #     # !!!! Until this part !!!!
    #
    #     return peptide

    @log_func
    def baseline_correction(self, input_file, random=True, threshold=5, i_baseline=0, include_negatives=False):
        '''This function takes the input_file DataFrame and substracts the baseline/noise channel from all other samples. The index of the
        baseline column is defaulted to 0. Set i_baseline=X to change baseline column.

        Threshold: After baseline substraction the remaining average signal has to be above threshold to be included. Parameter is set with threshold=X.
        This prevents very low remaining signal peptides to produce artificially high fold changes. Has to be determined empirically.

        First PSMs combined into peptides by 'Master Protein Accessions', 'Annotated Sequence', 'Modifications' and then baseline correction done, export
        peptide file
        '''

        random_float = np.random.RandomState(69)  # random seed for NaN, empty or 0 values.

        # channels and 'Modifications', 'Master Protein Accessions', 'Annotated Sequence' are required for further analysis for baseline correction.
        PSMs = input_file[self.channels + ['Modifications', 'Master Protein Accessions', 'Annotated Sequence']]

        # Group by 'Annotated Sequence' and 'Master Protein Accessions', 'Modifications' and aggregate the sum for each abundance column
        # the aim is to convert PSMs into peptide file.
        peptide = (
            PSMs.groupby(['Master Protein Accessions', 'Annotated Sequence', 'Modifications'])[self.channels]
                .agg('sum')
                .reset_index()
        )

        baseline_channel = self.channels[i_baseline]
        baseline = peptide[baseline_channel]
        peptide[self.channels] = peptide[self.channels].subtract(baseline, axis='index')
        peptide['Mean'] = peptide[self.channels].mean(axis=1)

        # !!!!! this part might change !!!!
        peptide = peptide.loc[peptide['Mean'] >= threshold]  # set S/N threshold for each PSM
        peptide = peptide.drop("Mean", axis=1)

        if (include_negatives == False and random == False):
            peptide[peptide < 0] = 0  # replace negative abundances with 0

        elif (include_negatives == False and random == True):
            for channel in self.channels:
                peptide[channel] = np.where(peptide[channel] < 0, random_float.random_sample(size=len(peptide)),
                                            peptide[channel])
        else:  # for other conditions we are not doing anything.
            pass

        # !!!! Until this part !!!!

        return peptide

    @log_func
    def protein_rollup(self, input_file, method='sum'):
        MPA=list([col for col in input_file.columns if 'Master Protein Accession' in col])[0]
        PSM_grouped=input_file.groupby(by=[MPA])
        result={}
        if method == 'sum':
            for group in PSM_grouped.groups:
                temp=PSM_grouped.get_group(group)
                sums=temp[self.channels].sum()
                result[group]=sums
        elif method == 'mean':
            for group in PSM_grouped.groups:
                temp=PSM_grouped.get_group(group)
                means = temp[self.channels].mean()
                result[group]=means

        elif method == 'median':
            for group in PSM_grouped.groups:
                temp=PSM_grouped.get_group(group)
                medians = temp[self.channels].median()
                result[group] = medians
                medians = temp[self.channels].median()
                result[group] = medians
        else:
            for group in PSM_grouped.groups:
                temp=PSM_grouped.get_group(group)
                medians = temp[self.channels].median()
                result[group] = medians
                sums=temp[self.channels].sum()
                result[group]=sums

        protein_df=pd.DataFrame.from_dict(result, orient='index',columns=self.channels)
        print("Protein rollup done!")
        return protein_df

    @log_func
    def log2(self, input):
        '''Modifies self.input_file and log2 transforms all TMT intensities.
        '''
        input[self.channels]=np.log2(input[self.channels])
        print("Normalization done")
        return input

    @log_func
    def return_file(self):
        return self.input_file

class plain_text_input:
    '''This class contains functions to analyze pSILAC data based on a plain text input file. The column names can be freely chosen, as
    long as all column names are unique. The different column identities are extracted by the column order:
    - Accession
    - Injection time (optional, is set in class init)
    - Modification
    - all following columns are assumed to contain TMT abundances
     '''
    def __init__(self, input, it_adj=True):
        '''Initialises class and extracts relevant columns.The different column identities are extracted by the column order:
        - Accession
        - Injection time (optional, set by it_adj parameter)
        - Modification
        - all following columns are assumed to contain TMT abundances 
         
        '''
        self.input_file = input
        self.input_columns = list(input.columns)
        if it_adj == True:
            self.abundances = self.input_columns[3:]
            self.mpa = self.input_columns[0]
            self.it_col = self.input_columns[1]
            self.modifications = self.input_columns[2]
        else:
            self.abundances = self.input_columns[2:]
            self.modifications = self.input_columns[1]
            self.mpa = self.input_columns[0]

    def get_channels(self, input):
        self.channels = [col for col in input.columns if 'Abundance:' in col]
        if self.channels == []:
            self.channels = [col for col in input.columns if 'Abundance' in col]

    def log_func(func):  # to show (log) what function is under usage!
        def wrapper(self, *args, **kwargs):
            print(f"Calling function: {func.__name__}")
            return func(self, *args, **kwargs)
        return wrapper

    @log_func
    def IT_adjustment(self, input):
        '''This function adjusts the input DataFrame stored in the class variable self.input_file for Ion injection times.
        '''
        self.get_channels(input)
        IT=self.it_col
        inject_times=input[IT]
        input[ self.channels]=input[ self.channels].divide(inject_times,axis=0)
        input[ self.channels]=input[ self.channels].multiply(1000)
        return input

    @log_func
    def filter_peptides(self, filtered_input):
        filtered_input = filtered_input[~filtered_input['Master Protein Accessions'].str.contains(';',na=False)]
        filtered_input = filtered_input[filtered_input['Contaminant'] == False]

        # this part removes the peptides with 0 signal/noise ratio (meaning that all the channels are 0)
        try:
            filtered_input = filtered_input.dropna(subset=['Average Reporter SN'])
            filtered_input = filtered_input[filtered_input['Average Reporter SN'] != 0]
        except:
            filtered_input = filtered_input.dropna(subset=['Average Reporter S/N'])
            filtered_input = filtered_input[filtered_input['Average Reporter S/N'] != 0]

        # this part removoes whole channel sum of 0, because we remove booster and rest might be 0 but Average Reporter SN is not 0
        # Calculate the row-wise sum for the specified columns
        row_sums = filtered_input[self.channels].sum(axis=1)
        # Filter out rows where the sum is zero
        filtered_input = filtered_input[row_sums != 0]

        # this part replace NaN values with 0
        filtered_input.fillna(0, inplace=True)

        return filtered_input

    @log_func
    def extract_heavy (self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all heavy labelled peptides. Naming of the 
        Modifications: Arg10: should contain Label, TMTK8, TMTproK8. Strings for modifications can be edited below for customisation.
        Returns heavy peptide DF
        '''
        modi=self.modifications
        Heavy_peptides=input[input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6|TMTproK4|TMTproK6',na=False)]
        print("Extraction Done","Extracted Heavy Peptides:", len(Heavy_peptides))
        return Heavy_peptides

    @log_func
    def extract_light (self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all light labelled peptides. Naming of the 
        Modifications: Arg10: should contain Label, TMTK8, TMTproK8. Strings for modifications can be edited below for customisation.

        Returns light peptide DF
        '''

        print("Extraction of light labelled peptides")
        modi=self.modifications

        light_peptides=input[~input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6|TMTproK4|TMTproK6',na=False)]

        print("Extraction Done","Extracted Heavy Peptides:", len(light_peptides))
        
        return light_peptides

    @log_func
    def baseline_correction(self,  input_file, random=True, threshold=5, i_baseline=0, include_negatives=False):
        '''This function takes the self.input_file DataFrame and substracts the baseline/noise channel from all other samples. The index of the
        baseline column is defaulted to 0. Set i_baseline=X to change baseline column.

        Threshold: After baseline substraction the remaining average signal has to be above threshold to be included. Parameter is set with threshold=X.
        This prevents very low remaining signal peptides to produce artificially high fold changes. Has to be determined empirically.

        Method: The method parameter sets the method for protein wollup quantification. Default is 'sum', which will sum all peptides for
        the corresponding protein. Alternatives are 'median' or 'mean'. If no or invalid input is given it uses 'sum'.

        Modifies self.input_file variable and returns a pandas df.
        '''

        random_float = np.random.RandomState(69)  # random seed for NaN, empty or 0 values.

        # channels and 'Modifications', 'Master Protein Accessions', 'Annotated Sequence' are required for further analysis for baseline correction.
        PSMs = input_file[self.channels + ['Modifications', 'Master Protein Accessions', 'Annotated Sequence']]

        # Group by 'Annotated Sequence' and 'Master Protein Accessions', 'Modifications' and aggregate the sum for each abundance column
        # the aim is to convert PSMs into peptide file.
        peptide = (
            PSMs.groupby(['Master Protein Accessions', 'Annotated Sequence', 'Modifications'])[self.channels]
                .agg('sum')
                .reset_index()
        )

        baseline_channel = self.channels[i_baseline]
        baseline = peptide[baseline_channel]
        peptide[self.channels] = peptide[self.channels].subtract(baseline, axis='index')
        peptide['Mean'] = peptide[self.channels].mean(axis=1)

        # !!!!! this part might change !!!!
        peptide = peptide.loc[peptide['Mean'] >= threshold]  # set S/N threshold for each PSM
        peptide = peptide.drop("Mean", axis=1)

        if (include_negatives == False and random == False):
            peptide[peptide < 0] = 0  # replace negative abundances with 0

        elif (include_negatives == False and random == True):
            for channel in self.channels:
                peptide[channel] = np.where(peptide[channel] < 0, random_float.random_sample(size=len(peptide)),
                                            peptide[channel])
        else:  # for other conditions we are not doing anything.
            pass

        # !!!! Until this part !!!!

        return peptide

    @log_func
    def TMM(self, input):
        '''This function implements TMM normalisation (Robinson & Oshlack, 2010, Genome Biology). It modifies the self.input_file class
        variable.
        '''
        channels=self.abundances
        input=input.dropna(subset=channels)
        input_trim=input[input[channels] < input[channels].quantile(.95)]
        print("Normalization")
        input_trim[channels]=input_trim[channels].divide(input_trim[channels[0]],axis=0)
        tm=np.argmin(trim_mean(input_trim[channels],0.25))
        summed=np.array(trim_mean(input_trim[channels], 0.25))
        minimum=summed[tm]
        norm_factors=summed/minimum
        input[channels]=input[channels].divide(norm_factors, axis=1)
        return input

    @log_func
    def total_intensity_normalisation(self, input):
        '''This function normalizes the self.input_file variable to the summed intensity of all TMT channels. It modifies the self.input_file
        to the updated DataFrame containing the normalized values.
        '''
        input_df = input.copy(deep=True)
        channels = [col for col in input_df.columns if 'Abundance:' in col]
        if not channels:
            channels = [col for col in input_df.columns if 'Abundance' in col]
        input_df = input_df.dropna(subset=channels)
        print("Normalization")
        minimum = np.argmin(input_df[channels].sum().values)
        summed = np.array(input_df[channels].sum().values)
        minimum = summed[minimum]
        norm_factors = summed / minimum
        input_df.loc[:, channels] = input_df[channels].divide(norm_factors, axis=1)
        print("Normalization done")
        return input_df

    @log_func
    def Median_normalisation(self, input):
        '''This function normalizes the self.input_file variable to the median of all individual TMT channels. It modifies the self.input_file
        to the updated DataFrame containing the normalized values.
        '''
        channels=self.abundances
        input=input.dropna(subset=channels)
        print("Normalization")
        minimum=np.argmin(input[channels].median().values)
        summed=np.array(input[channels].median().values)
        minimum=summed[minimum]
        norm_factors=summed/minimum
        input[channels]=input[channels].divide(norm_factors, axis=1)
        print("Normalization done")
        return input

    @log_func
    def protein_rollup(self, input_file, method='sum'):
        channels = self.abundances
        MPA = self.mpa
        PSM_grouped=input_file.groupby(by=[MPA])
        result={}
        if method == 'sum':
            for group in PSM_grouped.groups:
                temp=PSM_grouped.get_group(group)
                sums=temp[channels].sum()
                result[group]=sums
        elif method == 'mean':
            for group in PSM_grouped.groups:
                temp=PSM_grouped.get_group(group)
                means = temp[channels].mean()
                result[group]=means

        elif method == 'median':
            for group in PSM_grouped.groups:
                temp=PSM_grouped.get_group(group)
                medians = temp[channels].median()
                result[group] = medians
                medians = temp[channels].median()
                result[group] = medians
        else:
            for group in PSM_grouped.groups:
                temp=PSM_grouped.get_group(group)
                medians = temp[channels].median()
                result[group] = medians
                sums=temp[channels].sum()
                result[group]=sums

        protein_df=pd.DataFrame.from_dict(result, orient='index',columns=channels)

        return protein_df