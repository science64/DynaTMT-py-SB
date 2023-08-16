'''

    DynaTMT-py - a python package to process SILAC/TMT proteomics data
    Copyright (C) 2021  Kevin Klann

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
__version__ = "v2.2"
__maintainer__ = "Süleyman Bozkurt"
__email__ = "sbozkurt.mbg@gmail.com"
__date__ = '18.01.2021'
__update__ = '16.08.2023'

# from scipy import stats
from scipy.stats import trim_mean
import pandas as pd
import numpy as np
from numpy.random import random
import warnings
random_float = np.random.RandomState(69)
warnings.filterwarnings("ignore")

class PD_input:
    '''Class containing functions to analyze the default peptide/PSM output from ProteomeDiscoverer. All column names are assumed
    to be default PD output. If your column names do not match these assumed strings, you can modify them or use the plain_text_input
    class, that uses column order instead of names.
    '''
    def __init__(self):
        '''Initialises PD_input class with specified input file. The input file gets stored
        in the class variable self.input_file.
        '''
        # self.input_file = input # this has been changed, we are not using it anymore!
        pass

    def log_func(func):  # to show (log) what function is under usage!
        def wrapper(self, *args, **kwargs):
            print(f"Calling function: {func.__name__}")
            return func(self, *args, **kwargs)
        return wrapper

    @log_func
    def filter_peptides(self, filtered_input):
        filtered_input = filtered_input[~filtered_input['Master Protein Accessions'].str.contains(';',na=False)]
        filtered_input = filtered_input[filtered_input['Contaminant']==False]
        try:
            filtered_input = filtered_input.dropna(subset=['Average Reporter SN'])
            filtered_input = filtered_input[filtered_input['Average Reporter SN'] != 0]
        except:
            filtered_input = filtered_input.dropna(subset=['Average Reporter S/N'])
            filtered_input = filtered_input[filtered_input['Average Reporter S/N'] != 0]

        print("Filtering done!")
        return filtered_input

    @log_func
    def IT_adjustment(self, input):
        '''This function adjusts the input DataFrame stored in the class variable self.input_file for Ion injection times.
        Abundance channels should contain "Abundance:" string and injection time uses "Ion Inject Time" as used by ProteomeDiscoverer
        default output. For other column headers please refer to plain_text_input class.
        '''
        channels = [col for col in input.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in input.columns if 'Abundance' in col]
        IT=[col for col in input.columns if 'Ion Inject Time' in col]
        inject_times=input[IT[0]]
        input[channels]=input[channels].divide(inject_times,axis=0)
        input[channels]=input[channels].multiply(1000)
        print("IT adjustment done!")
        return input

    @log_func
    def extract_heavy (self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all heavy labelled peptides. Naming of the 
        Modifications: Arg10: should contain Label, TMTK8, TMTproK8. Strings for modifications can be edited below for customisation.
        Returns heavy peptide DF
        '''
        # print("Extraction of heavy labelled peptides")
        modi=list([col for col in input.columns if 'Modification' in col])
        modi=modi[0]
        '''Change Modification String here'''
        Heavy_peptides=input[input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6',na=False)]

        print("Extraction Done","Extracted Heavy Peptides:", len(Heavy_peptides))
        return Heavy_peptides

    @log_func
    def extract_light (self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all light labelled peptides. Naming of the 
        Modifications: Arg10: should contain Label, TMTK8, TMTproK8. Strings for modifications can be edited below for customisation.
        Returns light peptide DF
        '''
        # print("Extraction of light labelled peptides")
        modi=list([col for col in input.columns if 'Modification' in col])
        modi=modi[0]

        light_peptides=input[~input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6',na=False)]

        print("Extraction Done","Extracted Light Peptides:", len(light_peptides))
        return light_peptides

    @log_func
    def baseline_correction(self, input, threshold=5, i_baseline=0, method='sum'):
        '''This function takes the input_file DataFrame and substracts the baseline/noise channel from all other samples. The index of the
        baseline column is defaulted to 0. Set i_baseline=X to change baseline column.

        Threshold: After baseline substraction the remaining average signal has to be above threshold to be included. Parameter is set with threshold=X.
        This prevents very low remaining signal peptides to produce artificially high fold changes. Has to be determined empirically.

        Method: The method parameter sets the method for protein wollup quantification. Default is 'sum', which will sum all peptides for
        the corresponding protein. Alternatives are 'median' or 'mean'. If no or invalid input is given it uses 'sum'.

        Modifies self.input_file variable and returns a pandas df.
        '''
        # print("Baseline correction")
        channels=[col for col in input.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in input.columns if 'Abundance' in col]
        MPA = list([col for col in input.columns if 'Master Protein Accession' in col])
        MPA = MPA[0]
        protein_groups=input.groupby(by=[MPA],sort=False)
        results={}
        for group in protein_groups.groups:
            temp_data=protein_groups.get_group(group)
            temp_data = temp_data[channels]
            baseline_channel=channels[i_baseline]
            baseline=temp_data[baseline_channel]

            temp_data[channels]=temp_data[channels].subtract(baseline,axis='index')
            temp_data['Mean']=temp_data[channels].mean(axis=1)

            temp_data[temp_data < 0]=0 # replace negative abundances with 0
            temp_data=temp_data.loc[temp_data['Mean'] > threshold] # set S/N threshold for each PSM
            if method == 'sum':
                
                temp_data=temp_data[channels].sum()
                
            elif method == 'mean':
                temp_data=temp_data[channels].mean()
            elif method == 'median':
                temp_data=temp_data[channels].median()
            else:
                temp_data=temp_data[channels].sum()

            results[group]=temp_data

        print("Baseline correction done!")
        result_df=pd.DataFrame.from_dict(results, orient='index',columns=channels)
        self.input_file = result_df
        return result_df

    @log_func
    def baseline_correction_peptide_return(self, input_file, threshold=5, i_baseline=0, random=False, include_negatives=False):

        channels=[col for col in input_file.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in input_file.columns if 'Abundance' in col]
        MPA = [col for col in input_file.columns if 'Master Protein Accessions' in col]

        MPA = MPA[0]

        temp_data = input_file[channels]
        baseline_channel=channels[i_baseline]
        baseline=temp_data[baseline_channel]
        temp_data[channels]=temp_data[channels].subtract(baseline,axis='index')
        temp_data['Mean']=temp_data[channels].mean(axis=1)
        
        if (include_negatives==False and random==False):
            temp_data[temp_data < 0]=0 # replace negative abundances with 0
        elif (include_negatives == False and random ==True):
            for channel in channels:
                temp_data[channel] = np.where(temp_data[channel] < 0, random_float.random_sample(size=len(temp_data)),temp_data[channel])
        elif include_negatives == True:
            pass
        else:
            pass
        if threshold is None:
            pass
        else:
            temp_data=temp_data.loc[temp_data['Mean'] >= threshold] # set S/N threshold for each PSM
        input_file[channels]=temp_data[channels]
        return input_file

    @log_func
    def statistics(self, input):
        '''This function provides summary statistics for quality control assessment from Proteome Discoverer Output.
        '''
        numeric_without_tmt_channels=[col for col in input.columns if not 'Abundance:' in col]
        if numeric_without_tmt_channels == []:
            channels = [col for col in input.columns if 'Abundance' in col]
        print(input[numeric_without_tmt_channels].describe(include=[np.number]))
        return(input[numeric_without_tmt_channels].describe(include=[np.number]))

    @log_func
    def TMM(self, input):
        '''This function implements TMM normalisation (Robinson & Oshlack, 2010, Genome Biology). It modifies the self.input_file class
        variable.
        '''
        channels=[col for col in input.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in input.columns if 'Abundance' in col]
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
    def chunks(self,l, n):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i + n]

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
        minimum = np.argmin(input_df[channels].sum().values)
        summed = np.array(input_df[channels].sum().values)
        minimum = summed[minimum]
        norm_factors = summed / minimum
        input_df.loc[:, channels] = input_df[channels].divide(norm_factors, axis=1)
        print("Total intensity normalisation done!")
        return input_df

    @log_func
    def Median_normalisation(self, input):
        '''This function normalizes the self.input_file variable to the median of all individual TMT channels. It modifies the self.input_file
        to the updated DataFrame containing the normalized values.
        '''
        channels=[col for col in input.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in input.columns if 'Abundance' in col]
        input=input.dropna(subset=channels)
        minimum=np.argmin(input[channels].median().values)
        summed=np.array(input[channels].median().values)
        minimum=summed[minimum]
        norm_factors=summed/minimum
        input[channels]=input[channels].divide(norm_factors, axis=1)
        print("Median normalisation done!")
        return input

    @log_func
    def protein_rollup(self, input_file, method='sum'):
        channels = [col for col in input_file.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in input_file.columns if 'Abundance' in col]
        MPA=list([col for col in input_file.columns if 'Master Protein Accession' in col])
        MPA=MPA[0]
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
        print("Protein rollup done!")
        return protein_df

    @log_func
    def sum_peptides_for_proteins(self, input_file):
        '''This function takes a peptide/PSM level DataFrame stored in self.input_file and performs Protein quantification rollup based
        on the sum of all corresponding peptides.

        Returns a Protein level DataFrame and modifies self.input_file
        '''
        print('Calculate Protein quantifications from PSM!')
        
        channels = [col for col in input_file.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in input_file.columns if 'Abundance' in col]
        MPA=list([col for col in input_file.columns if 'Master Protein Accession' in col])
        MPA=MPA[0]
        PSM_grouped=input_file.groupby(by=[MPA])
        result={}
        for group in PSM_grouped.groups:
            temp=PSM_grouped.get_group(group)
            sums=temp[channels].sum()
            result[group]=sums
        
        protein_df=pd.DataFrame.from_dict(result, orient='index',columns=channels)
        print("Protein rollup done!")
        return protein_df[channels]

    @log_func
    def log2(self, input):
        '''Modifies self.input_file and log2 transforms all TMT intensities.
        '''
        channels=[col for col in input.columns if 'Abundance:' in col]
        if channels == []:
            channels = [col for col in input.columns if 'Abundance' in col]
        input[channels]=np.log2(input[channels])
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

    def log_func(func):  # to show (log) what function is under usage!
        def wrapper(self, *args, **kwargs):
            print(f"Calling function: {func.__name__}")
            return func(self, *args, **kwargs)
        return wrapper

    @log_func
    def IT_adjustment(self, input):
        '''This function adjusts the input DataFrame stored in the class variable self.input_file for Ion injection times.
        
        '''
        print("IT adjustment")
        channels = self.abundances
        IT=self.it_col
        inject_times=input[IT]
        input[channels]=input[channels].divide(inject_times,axis=0)
        input[channels]=input[channels].multiply(1000)
        print("Done")
        return input

    @log_func
    def extract_heavy (self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all heavy labelled peptides. Naming of the 
        Modifications: Arg10: should contain Label, TMTK8, TMTproK8. Strings for modifications can be edited below for customisation.
        Returns heavy peptide DF
        '''

        print("Extraction of heavy labelled peptides")
        modi=self.modifications

        Heavy_peptides=input[input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6',na=False)]

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

        light_peptides=input[~input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6',na=False)]

        print("Extraction Done","Extracted Heavy Peptides:", len(light_peptides))
        
        return light_peptides

    @log_func
    def baseline_correction(self, input, threshold=5, i_baseline=0, method='sum'):
        '''This function takes the self.input_file DataFrame and substracts the baseline/noise channel from all other samples. The index of the
        baseline column is defaulted to 0. Set i_baseline=X to change baseline column.

        Threshold: After baseline substraction the remaining average signal has to be above threshold to be included. Parameter is set with threshold=X.
        This prevents very low remaining signal peptides to produce artificially high fold changes. Has to be determined empirically.

        Method: The method parameter sets the method for protein wollup quantification. Default is 'sum', which will sum all peptides for
        the corresponding protein. Alternatives are 'median' or 'mean'. If no or invalid input is given it uses 'sum'.

        Modifies self.input_file variable and returns a pandas df.
        '''
        
        print("Baseline correction")
        channels=self.abundances
        MPA = self.mpa
        
        protein_groups=input.groupby(by=[MPA],sort=False)
        results={}
        for group in protein_groups.groups:
            temp_data=protein_groups.get_group(group)
            temp_data = temp_data[channels]
            baseline_channel=channels[i_baseline]
            baseline=temp_data[baseline_channel]

            temp_data[channels]=temp_data[channels].subtract(baseline,axis='index')
            temp_data['Mean']=temp_data[channels].mean(axis=1)

            temp_data[temp_data < 0]=0 # replace negative abundances with 0
            temp_data=temp_data.loc[temp_data['Mean'] > threshold] # set S/N threshold for each PSM
            if method == 'sum':
                temp_data=temp_data[channels].sum()
            elif method == 'mean':
                temp_data=temp_data[channels].mean()
            elif method == 'median':
                temp_data=temp_data[channels].median()
            else:
                temp_data=temp_data[channels].sum()

            results[group]=temp_data
        print("Baseline correction done")
        result_df=pd.DataFrame.from_dict(results, orient='index',columns=channels)
        self.input_file = result_df
        return result_df

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

    @log_func
    def sum_peptides_for_proteins(self, input):
        '''This function takes a peptide/PSM level DataFrame stored in self.input_file and performs Protein quantification rollup based
        on the sum of all corresponding peptides.

        Returns a Protein level DataFrame and modifies self.input_file
        '''
        print('Calculate Protein quantifications from PSM')
        channels = self.abundances
        MPA=self.mpa
        PSM_grouped=input.groupby(by=[MPA])
        result={}
        for group in PSM_grouped.groups:
            temp=PSM_grouped.get_group(group)
            sums=temp[channels].sum()
            result[group]=sums
        
        protein_df=pd.DataFrame.from_dict(result, orient='index',columns=channels)
        
        print("Combination done")
        return protein_df[channels]

    @log_func
    def baseline_correction_peptide_return(self, input_file, threshold=5, i_baseline=0, random=False, include_negatives=False):
        print("Baseline correction")
        channels=self.abundances
        MPA = self.mpa
        protein_groups=input_file.groupby(by=[MPA],sort=False)
        results={}

        temp_data = input_file[channels]
        baseline_channel=channels[i_baseline]
        baseline=temp_data[baseline_channel]
        temp_data[channels]=temp_data[channels].subtract(baseline,axis='index')
        temp_data['Mean']=temp_data[channels].mean(axis=1)
        
        if (include_negatives==False and random==False):
            temp_data[temp_data < 0]=0 # replace negative abundances with 0
        elif (include_negatives == False and random ==True):
            for channel in channels:
                temp_data[channel] = np.where(temp_data[channel] < 0, random_float.random_sample(size=len(temp_data)),temp_data[channel])
        elif include_negatives == True:
            pass
        else:
            pass
        temp_data=temp_data.loc[temp_data['Mean'] > threshold] # set S/N threshold for each PSM
        input_file[channels]=temp_data[channels]
        return input_file

    @log_func
    def return_file(self):
        return self.input_file