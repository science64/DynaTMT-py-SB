'''

    DynaTMT-py - a python package to process SILAC/TMT proteomics data
    Copyright (C) 2021  Kevin Klann - 2024 Süleyman Bozkurt

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
__version__ = "v2.9.2"
__maintainer__ = "Süleyman Bozkurt"
__email__ = "sbozkurt.mbg@gmail.com"
__date__ = '18.01.2021'
__update__ = '03.06.2024'

from scipy.stats import trim_mean
import pandas as pd
import numpy as np
from numpy.random import random
import warnings
import re
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
        # return unnormalized abundances column names (channels)
        channels = [col for col in input.columns
                    if 'abundance' in col.lower() and 'normaliz' not in col.lower()]
        return channels

    def log_func(func):  # to show (log) what function is under usage!
        def wrapper(self, *args, **kwargs):
            print(f"Calling function: {func.__name__}")
            return func(self, *args, **kwargs)
        return wrapper

    @log_func
    def filter_peptides(self, filtered_input):

        # remove shared peptides (peptides with more than one protein accession)
        filtered_input = filtered_input[filtered_input['Quan Info'] != 'NotUnique']
        # remove contaminants
        filtered_input = filtered_input[filtered_input['Contaminant'] == False]

        # Remove rows where any of the specified 'channels' columns have at least one NA value
        filtered_input = filtered_input.dropna(subset=self.channels)

        return filtered_input

    @log_func
    def filter_PSMs(self, filtered_input):

        ### Apply this filter to the PSMs only ###
        # the idea is to find isolation interference column and filter out the PSMs with isolation interference > 50%
        isolation_interference_col = None
        # Loop through the columns and find the one that contains 'isolation interference'
        for col in filtered_input.columns:
            if 'isolation interference' in col.lower():
                isolation_interference_col = col
                break  # Stop the loop once the column is found

        # Apply the filter using the identified column name
        if isolation_interference_col:
            filtered_input = filtered_input[
                filtered_input[isolation_interference_col] < 50]

        # remove shared peptides (peptides with more than one protein accession)
        filtered_input = filtered_input[~filtered_input['Master Protein Accessions'].str.contains(';', na = False)]

        # remove empty accessions (if any)
        filtered_input = filtered_input.dropna(subset=['Master Protein Accessions'])

        # remove contaminants
        filtered_input = filtered_input[filtered_input['Contaminant'] == False]

        # Remove rows where any of the specified 'channels' columns have at least one NA value
        filtered_input = filtered_input.dropna(subset=self.channels)

        return filtered_input

    @log_func
    def IT_adjustment(self, input):
        '''
        This function is only important for MS2 measuments if needed. (OPTIONAL)
        This function adjusts the input DataFrame stored in the class variable self.input_file for Ion injection times.
        Abundance channels should contain "Abundance:" string and injection time uses "Ion Inject Time" as used by ProteomeDiscoverer
        default output. For other column headers please refer to plain_text_input class.
        '''
        IT=[col for col in input.columns if 'Ion Inject Time' in col]
        inject_times=input[IT[0]]
        input[self.channels]=input[self.channels].divide(inject_times, axis=0)
        input[self.channels]=input[self.channels].multiply(1000)
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
    def extract_heavy(self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all heavy labelled peptides. Naming of the 
            Modifications: Heavy Arginine should contain: Label
                           Heavy Lysine should contain: TMTK4, TMTK6, TMTK8, TMTproK4, TMTproK6, TMTproK8.
                           Strings for modifications can be edited below for customisation.
            Returns heavy peptide DF
        '''
        modi = list([col for col in input.columns if 'Modification' in col])[0]
        '''Change Modification String here'''
        Heavy_peptides = input[input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6|TMTproK4|TMTproK6', na=False)]
        print("Extraction Done","Extracted Heavy PSMs/Peptides:", len(Heavy_peptides))
        return Heavy_peptides

    @log_func
    def extract_light(self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all light labelled peptides. Naming of the 
            Modifications: Light Arginine should NOT contain: Label
                           Light Lysine should NOT contain: TMTK4, TMTK6, TMTK8, TMTproK4, TMTproK6, TMTproK8.
                           Strings for modifications can be edited below for customisation.
            Returns light peptide DF
        '''
        modi = list([col for col in input.columns if 'Modification' in col])[0]
        light_peptides = input[~input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6|TMTproK4|TMTproK6',na=False)]
        print("Extraction Done","Extracted Light PSMs/Peptides:", len(light_peptides))
        return light_peptides

    @log_func
    def PSMs_to_Peptide(self, input):
        '''
        - This function takes PSMs as the input file and groups the PSMs into peptides by the
        Annotated sequence, Modifications and Master Protein Accessions. It then aggregates the TMT channels.
        Finally, it returns a DataFrame containing the peptide sequences, modifications, accession ID
        and the aggregated TMT channels from peptides.

        - If the required columns are not found, the function will try to find the column that contains the theoretical MH+ values.
        '''

        try:
            peptides = input.groupby(
                ['Annotated Sequence', 'Modifications', 'Master Protein Accessions'])[
                self.channels].sum().reset_index()
        except:
            # if the dataframe does not contain the required columns, we need to find the column that contains the theoretical MH+ values
            '''
            First we need to identify the column name that contains the theoretical MH+ values. 
            This is necessary because PSMs matches and combining them into peptides.
            '''
            # List of exact column names to search for
            exact_col_names = ['Theo. MH+ [Da]', 'Theo MHplus in Da']

            # Variable to hold the found column name or None if not found
            found_col_name = None

            # First, try to find an exact match from the list
            for col in exact_col_names:
                if col in input.columns:
                    found_col_name = col
                    break

            # If no exact match is found, use regex to find a matching column
            if found_col_name is None:
                standard_col_name = 'Theo_MHplus_Da'  # This will be the new standard column name
                for col in input.columns:
                    if re.search(r'theo.? mh\+? \[?da\]?', col, re.IGNORECASE):
                        input.rename(columns={col: standard_col_name}, inplace=True)
                        found_col_name = standard_col_name
                        break

            # Now, use the found column name in your groupby if it's not None
            if found_col_name is not None:
                peptides = input.groupby(['Master Protein Accessions', found_col_name])[self.channels].sum().reset_index()
            else:
                print("The required Theo MH+ column was not found in the DataFrame.")
                return None

        return peptides

    @log_func
    def baseline_correction(self, input_file, threshold=5, i_baseline=0, random=True):  # include_negatives=False, Because there is no use of it!
        """
        This function takes the input_file DataFrame and substracts the baseline/noise channel from all other samples. The index of the
        baseline column is defaulted to 0. Set i_baseline=X to change baseline column.

        if random is True, it will replace the negative values with random values between 0 and 1.

        Threshold: After baseline substraction the remaining average signal has to be above threshold to be included. Parameter is set with threshold=X.
        This prevents very low remaining signal peptides to produce artificially high fold changes. Has to be determined empirically.

        It can identify the file is PSMs or Peptides by the column name. Then it will do the baseline correction for PSMs or Peptides.
        It will convert PSMs into Peptides by sum all the same ('Master Protein Accessions', 'Annotated Sequence', 'Modifications') at the last step.
        """

        # determine input file is peptide or PSMs
        if 'PSMs Peptide ID' in input_file.columns:
            decision = 'PSMs'
        elif any('Peptide Group ID' in column for column in input_file.columns):
            decision = 'Peptides'
        else:
            decision = 'Unknown'

        print('[#] Decision of this file is:', decision)

        baseline_channel = self.channels[i_baseline]
        baseline = input_file[baseline_channel]
        input_file[self.channels] = input_file[self.channels].subtract(baseline, axis='index')

        if random == True:
            random_float = np.random.RandomState(69)
            # Step 1: Replace all negative values with 0 in all channels
            for channel in self.channels:
                input_file[channel][input_file[channel] < 0] = 0

            # then again I need to sum all the values across the channels and remove the rows where the sum is 0
            # Step 2: Calculate the sum of values across specified channels
            input_file['sum_abundances'] = input_file[self.channels].sum(axis=1)

            # Step 3: Filter out rows where 'sum_abundances' is zero
            input_file = input_file[input_file['sum_abundances'] != 0]

            # Step 4: Drop the 'sum_abundances' column
            input_file.drop(columns=['sum_abundances'], inplace=True)

            # step 5: Remove baseline channel from the channels list
            channels_final = self.channels.copy()
            channels_final.remove(baseline_channel)

            # Step 6: Assign random values between 0 and 1 to the 0 values and round to two decimal places
            for channel in channels_final:
                input_file[channel] = np.where(input_file[channel] <= 0,
                                            random_float.random_sample(size=len(input_file[channel])),
                                            input_file[channel])

                # # Ensure all values in the copy are rounded to two decimal places
                # input_file[channel] = np.round(input_file[channel], 2)

        else:
            # all the negative values will be replaced by 0 if not random
            for channel in self.channels:
                input_file[channel][input_file[channel] < 0] = 0

        # Step 1: Calculate the mean across specified channels
        input_file['Mean'] = input_file[self.channels].mean(axis=1)
        # Step 2: Filter the DataFrame based on the 'Mean' column with the given threshold
        input_file = input_file.loc[input_file['Mean'] >= threshold]  # set S/N threshold for each PSM
        # Step 3: Drop the 'Mean' column
        input_file = input_file.drop("Mean", axis=1)

        # print('PSMs_data_2 with %s threshold: %s rows x %s columns' % (
        # threshold, input_file.shape[0], input_file.shape[1]))

        if decision == 'PSMs':
            # Group by 'Annotated Sequence' and 'Master Protein Accessions', 'Modifications' and aggregate the sum for each abundance column
            # the aim is to convert PSMs into peptide file.

            peptides = self.PSMs_to_Peptide(input_file)  # convert PSMs into Peptides
        else:
            peptides = input_file.copy()

        return peptides

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
        return protein_df

    @log_func
    def log2(self, input):
        '''Modifies self.input_file and log2 transforms all TMT intensities.
        '''
        input[self.channels]=np.log2(input[self.channels])
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
    def extract_heavy(self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all heavy labelled peptides. Naming of the 
        Modifications: Arg10: should contain Label, TMTK8, TMTproK8. Strings for modifications can be edited below for customisation.
        Returns heavy peptide DF
        '''
        modi=self.modifications
        Heavy_peptides=input[input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6|TMTproK4|TMTproK6',na=False)]
        print("Extraction Done","Extracted Heavy Peptides:", len(Heavy_peptides))
        return Heavy_peptides

    @log_func
    def extract_light(self, input):
        '''This function takes the class variable self.input_file dataframe and extracts all light labelled peptides. Naming of the 
        Modifications: Arg10: should contain Label, TMTK8, TMTproK8. Strings for modifications can be edited below for customisation.

        Returns light peptide DF
        '''

        print("Extraction of light labelled peptides")
        modi=self.modifications

        light_peptides=input[~input[modi].str.contains('TMTK8|Label|TMTproK8|TMTK4|TMTK6|TMTproK4|TMTproK6',na=False)]

        print("Extraction Done","Extracted Heavy Peptides:", len(light_peptides))
        
        return light_peptides

    def PSMs_to_Peptide(self, input):
        '''This function takes PSMs as the input file and groups the PSMs into peptides by the
        Annotated sequence, Modifications and Master Protein Accessions. It then aggregates the TMT channels.
        It then finally returns a DataFrame containing the peptide sequences, modifications, accession ID
        and the aggregated TMT channels.
        '''

        try:
            peptides = input.groupby(
                ['Annotated Sequence', 'Modifications', 'Master Protein Accessions'])[
                self.channels].sum().reset_index()
        except:
            # if the dataframe does not contain the required columns, we need to find the column that contains the theoretical MH+ values
            '''
            First we need to identify the column name that contains the theoretical MH+ values. 
            This is necessary because PSMs matches and combining them into peptides.
            '''
            # List of exact column names to search for
            exact_col_names = ['Theo. MH+ [Da]', 'Theo MHplus in Da']

            # Variable to hold the found column name or None if not found
            found_col_name = None

            # First, try to find an exact match from the list
            for col in exact_col_names:
                if col in input.columns:
                    found_col_name = col
                    break

            # If no exact match is found, use regex to find a matching column
            if found_col_name is None:
                standard_col_name = 'Theo_MHplus_Da'  # This will be the new standard column name
                for col in input.columns:
                    if re.search(r'theo.? mh\+? \[?da\]?', col, re.IGNORECASE):
                        input.rename(columns={col: standard_col_name}, inplace=True)
                        found_col_name = standard_col_name
                        break

            # Now, use the found column name in your groupby if it's not None
            if found_col_name is not None:
                peptides = input.groupby(['Master Protein Accessions', found_col_name])[self.channels].sum().reset_index()
            else:
                print("The required Theo MH+ column was not found in the DataFrame.")
                return None

        return peptides

    @log_func
    def baseline_correction(self, input_file, threshold=5, i_baseline=0,
                            random=True):  # include_negatives=False, Because there is no use of it!
        '''This function takes the input_file DataFrame and substracts the baseline/noise channel from all other samples. The index of the
        baseline column is defaulted to 0. Set i_baseline=X to change baseline column.

        if random is True, it will replace the negative values with random values between 0 and 1.

        Threshold: After baseline substraction the remaining average signal has to be above threshold to be included. Parameter is set with threshold=X.
        This prevents very low remaining signal peptides to produce artificially high fold changes. Has to be determined empirically.

        It can identify the file is PSMs or Peptides by the column name. Then it will do the baseline correction for PSMs or Peptides.
        It will convert PSMs into Peptides by sum all the same ('Master Protein Accessions', 'Annotated Sequence', 'Modifications') at the last step.
        '''

        # determine input file is peptide or PSMs
        if 'PSMs Peptide ID' in input_file.columns:
            decision = 'PSMs'
        elif any('Peptide Group ID' in column for column in input_file.columns):
            decision = 'Peptides'
        else:
            decision = 'Unknown'

        print('[#] Decision of this file is:', decision)

        baseline_channel = self.channels[i_baseline]
        baseline = input_file[baseline_channel]
        input_file[self.channels] = input_file[self.channels].subtract(baseline, axis='index')

        if random == True:
            random_float = np.random.RandomState(69)
            # Step 1: Replace all negative values with 0 in all channels
            for channel in self.channels:
                input_file[channel][input_file[channel] < 0] = 0

            # then again I need to sum all the values across the channels and remove the rows where the sum is 0
            # Step 2: Calculate the sum of values across specified channels
            input_file['sum_abundances'] = input_file[self.channels].sum(axis=1)

            # Step 3: Filter out rows where 'sum_abundances' is zero
            input_file = input_file[input_file['sum_abundances'] != 0]

            # Step 4: Drop the 'sum_abundances' column
            input_file.drop(columns=['sum_abundances'], inplace=True)

            # step 5: Remove baseline channel from the channels list
            channels_final = self.channels.copy()
            channels_final.remove(baseline_channel)

            # Step 5: Assign random values between 0 and 1 to the 0 values and round to two decimal places
            for channel in channels_final:
                input_file[channel] = np.where(input_file[channel] <= 0,
                                               random_float.random_sample(size=len(input_file[channel])),
                                               input_file[channel])

                # Ensure all values in the copy are rounded to two decimal places
                input_file[channel] = np.round(input_file[channel], 2)

        else:
            # all the negative values will be replaced by 0 if not random
            for channel in self.channels:
                input_file[channel][input_file[channel] < 0] = 0

        # Step 1: Calculate the mean across specified channels
        input_file['Mean'] = input_file[self.channels].mean(axis=1)
        # Step 2: Filter the DataFrame based on the 'Mean' column with the given threshold
        input_file = input_file.loc[input_file['Mean'] >= threshold]  # set S/N threshold for each PSM
        # Step 3: Drop the 'Mean' column
        input_file = input_file.drop("Mean", axis=1)

        if decision == 'PSMs':
            # Group by 'Annotated Sequence' and 'Master Protein Accessions', 'Modifications' and aggregate the sum for each abundance column
            # the aim is to convert PSMs into peptide file.
            peptides = self.PSMs_to_Peptide(input_file)  # convert PSMs into Peptides
        else:
            peptides = input_file.copy()

        return peptides

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