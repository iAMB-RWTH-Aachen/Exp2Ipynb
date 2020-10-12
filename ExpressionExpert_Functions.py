"""
Functions for the ExpressionExpertIpynb workflow.

Author: Ulf Liebal
Contact: ulf.liebal@rwth-aachen.de
Date: 2020, May
"""
###########################################################################
###########################################################################
#
# Data preparation and system set-up
#
###########################################################################
###########################################################################

def init_Exp2(config_File):
    '''
    Parameter initialization for ExpressionExpertIpynb.
    
    Input:
        config_File:    string; file address of the configuration file with parameters for statistical analysis and regression.
        
    Output:
        Name_Dict:      dictionary; parameter dictionary
    '''
    Name_Dict = dict()
    with open(config_File) as Conf:
        myline = Conf.read().splitlines()
        for line in myline:
            if not line.startswith('#'):
                (key, val) = line.split(':', 1)
                Name_Dict[str(key.strip())] = val.strip()
#     Data_File = Name_Dict['Data_File']
#     # extract the filename for naming of newly generated files
#     File_Base = Name_Dict['File_Base']
#     # the generated files will be stored in a subfolder with custom name
#     Data_Folder = 'data-{}'.format(File_Base)    # column name of expression values
#     Y_Col_Name = eval(Name_Dict['Y_Col_Name'])
#     # Extracting entropy cutoff for removal of non-informative positions
#     Entropy_cutoff = float(Name_Dict['Entropy_cutoff'])
#     # figure file type
#     Fig_Type = Name_Dict['Figure_Type']
#     # Figure font size
#     FigFontSize = Name_Dict['Figure_Font_Size']
    make_DataDir(Name_Dict)
    
    return Name_Dict

def make_DataDir(Name_Dict):
    '''
    Set-up of directory for data storage.
    
    Input:
            Name_Dict:      dictionary; contains parameters as defined in the configuration file 'congig.txt'

    Output: 
            None
    '''
    import os
    
    Data_File = Name_Dict['Data_File']
    # extract the filename for naming of newly generated files
    File_Base = Data_File.split('.')[0]
    # the generated files will be stored in a subfolder with custom name
    Data_Folder = 'data-{}'.format(File_Base)
    try:
        os.mkdir(Data_Folder)
        print('Data directory ', Data_Folder, 'created.')
    except FileExistsError:
        print('Already existent data directory ', Data_Folder, '.')

###########################################################################
###########################################################################

def Data_Src_Load(Name_Dict):
    '''
    Data loading of Excel sheet to pandas data frame. The dictionary defines the file name, and columns for sequence and expression strength. The sequence is converted into a one-hot-, label encoding (labels 0-3) and the GC-content is added. 
    
    Input:
            Name_Dict:      dictionary; contains parameters as defined in the configuration file 'congig.txt'
 
    Output:
            SeqDat:         dataframe; sequence in letter, label, and one-hot format, GC-content, expression strength and sequence IDs
    '''
#     import os
    import numpy as np
    import pandas as pd
    from sklearn.preprocessing import StandardScaler
#     from ExpressionExpert_Functions import list_integer, list_onehot

    DataPath = Name_Dict['Data_File']
    Seq_Col = Name_Dict['Sequence_column']
    ID_Col_Name = Name_Dict['ID_Col_Name']
    Y_Col_Name = eval(Name_Dict['Y_Col_Name'])
    ML_Regressor = Name_Dict['ML_Regressor']
    Response_Value = eval(Name_Dict['Response_Value'])
    
    if 'Remove_Outlier' in Name_Dict:
        Remove_Outlier = eval(Name_Dict['Remove_Outlier'])
    else:
        Remove_Outlier = False

    
    SeqDat = pd.read_csv(DataPath, delimiter=',|;', engine='python')
    # decision whether to delete samples identified as outlier
    # outlier are defined as farther than 3x of standard deviation from the mean value
    if Remove_Outlier:
        myOutlier=list()
        for Measurement in Y_Col_Name:
            myOutlier.append(DetectOutlier(SeqDat[Measurement]))
        myOutlier = np.array(myOutlier).ravel()
        print('Following outliers were detected: ID: {}, Value: {}'.format(SeqDat.loc[myOutlier, ID_Col_Name].values,SeqDat.loc[myOutlier, Y_Col_Name].values))
        SeqDat.drop(myOutlier, inplace=True)
        SeqDat.reset_index(inplace=True)
    
    # decision whether additional data is generated from statistics, this requires mean standard deviation and sample size
    Stats2Samples = eval(Name_Dict['Stats2Samples'])
    if Stats2Samples:
        SeqDat = SeqDat_StatsExpand(SeqDat, Name_Dict)
    
    # testing whether all sequences have the same length
    # sequences deviating from the median length are deleted
    SLvec = list()
    for item in SeqDat[Seq_Col]:
        SLvec.append(len(item))
    SequenceLengthError = np.where(np.array(SLvec)!=np.median(SLvec))
    SeqDat.drop(SeqDat.index[SequenceLengthError], inplace=True)
#     SeqDat.reset_index(inplace=True)
    
    SeqDat['Sequence_label-encrypted'] = list_integer(SeqDat[Seq_Col].str.upper())
    SeqDat['Sequence_letter-encrypted'] = SeqDat[Seq_Col].str.upper()
    SeqDat['Sequence'] = list_onehot(SeqDat['Sequence_label-encrypted'])
    # GC content calculation
    NuclSum = [np.sum(SeqDat['Sequence'][i], axis=0) for i in range(len(SeqDat))]
    GCcont = np.sum(np.delete(np.vstack(NuclSum),[0,3],1), axis=1)/np.sum(np.vstack(NuclSum)[0])
    SeqDat['GC-content'] = GCcont
        
    # Categorization of expression
    # Extracting the last letter of the abbreviation of the ML approach.
    # If the last letter is a 'R', a regression is chosen, if 'C' a classification is chosen.
    # If classification is used, the expression results are categorized
    if Response_Value == 0:
        # using standardized data with zero mean and unit variance
        print('Standardizing of expression.')
        for Measurement in Y_Col_Name:
            myScaler = StandardScaler()
            TargetName = '{}_ML'.format(Measurement)
            myData = SeqDat[Measurement].values.reshape(-1, 1)
            SeqDat[TargetName] = myScaler.fit_transform(myData)
        
    elif Response_Value == 1:
        print('Original values of expression.')
        for Measurement in Y_Col_Name:
            TargetName = '{}_ML'.format(Measurement)
            SeqDat[TargetName] = SeqDat[Measurement]
            
    elif Response_Value > 1:
        print('Categorization of expression.')
        BinNum = Response_Value
        for Measurement in Y_Col_Name:
            CatName = '{}_ML'.format(Measurement)
            SeqDat[CatName], bins = pd.qcut(SeqDat[Measurement], q=BinNum, labels=range(BinNum), retbins=True)
        SeqDat.drop_duplicates('Sequence_letter-encrypted', inplace=True)   
        print('The expression values were sorted into the following bins: {}'.format(bins))
        
    else:
        print('Response value parameter must be an positive integer number')
    SeqDat.reset_index(inplace=True)
    

    return SeqDat

###########################################################################
###########################################################################

###########################################################################
# The following functions serve for the one-hot encoding
# It is derived from:
#    https://machinelearningmastery.com/how-to-one-hot-encode-sequence-data-in-python/
def list_integer(SeqList):
    '''
    Integer label encryption of base letter sequence.
    
    Input:
        SeqList:       list; letter sequence of nucleotide bases A,C,G,T
        
    Output:
        IntegerList:   list; numerical label encryption of letters A-0, C-1, G-2, T-3
    '''
    from sklearn.preprocessing import LabelEncoder
    
    alphabet = 'ACGT'
    LabEnc = LabelEncoder().fit(list(alphabet))
    IntegerList = [LabEnc.transform(list(i)) for i in SeqList]
    return IntegerList    
            
def list_onehot(IntegerList):
    '''
    Generate one-hot encoding.
    
    Input:
        IntegerList:    list, label encrypted sequence
        
    Output:
        OneHotList:     list, one-host encoded sequence
    '''
    OneHotList = list()
    for integer_encoded in IntegerList:    
        onehot_encoded = list()
        for value in integer_encoded:
            letter = [0 for _ in range(4)]
            letter[value] = 1
            onehot_encoded.append(letter)
        OneHotList.append(onehot_encoded)
    return OneHotList

###########################################################################
###########################################################################

def ExpressionScaler(SeqDat, Name_Dict):
    '''
    Scaling of expression values to zero mean and unit variance.
    
    Input:
        SeqDat:         dataframe; sequence in letter, label, and one-hot format, GC-content, expression strength and sequence IDs
        Name_Dict:      dictionary; contains parameters as defined in the configuration file 'congig.txt'
        
     Output:
         SeqDat:      dataframe; like input with added scaled measurements
         Expr_Scaler: dictionary; contains the scaler functions
    '''
    from sklearn.preprocessing import StandardScaler
    import copy
    
    Y_Col_Name = eval(Name_Dict['Y_Col_Name'])
    Measurement_Number = int(Name_Dict['Library_Expression'])
    
    Expr_Scaler = dict()
    Expr_Scaler_n = StandardScaler() 
    for idx in range(Measurement_Number):
        Column_Name = '{}_ML'.format(Y_Col_Name[idx])
        myData = SeqDat[Y_Col_Name[idx]].values.reshape(-1, 1)
        SeqDat[Column_Name] = Expr_Scaler_n.fit_transform(myData)
        Scaler_Name = '{}_Scaler'.format(Y_Col_Name[idx])
        Expr_Scaler[Scaler_Name] = copy.deepcopy(Expr_Scaler_n)
     
    return SeqDat, Expr_Scaler
    
###########################################################################
###########################################################################
    
def Insert_row_(row_number, df, row_value): 
    '''
    The function allows the insertion of rows in a dataframe. The index counting 
    is reversed, the last element is set to -1 and each element to the top is 
    decreased by 1.
    
    Input:
        row_number: numpy integer/vector; index from top to insert line
        df:         dataframe
        row_value:  numpy vector/matrix; contains numbers being inserted
    '''
    import pandas as pd
    
    for idx, element in zip(row_number, row_value):
        # Slice the upper half of the dataframe 
        df1 = df[0:idx] 

        # Store the result of lower half of the dataframe 
        df2 = df[idx:] 

        # Inser the row in the upper half dataframe 
        df1.loc[idx]=element 

        # Concat the two dataframes 
        df = pd.concat([df1, df2]) 

        # Reassign the index labels 
        df.index = [*range(df.shape[0])] 

    # Return the updated dataframe 
    df.index = [*range(-df.shape[0],0,1)] 
    return df

###########################################################################
###########################################################################

def SeqDat_StatsExpand(SeqDat, Name_Dict):
    '''
    Generates separate samples from statistical information including mean, standard deviation and sample size.

        SeqDat:         dataframe; sequence in letter, label, and one-hot format, GC-content, expression strength and sequence IDs
        Name_Dict:      dictionary; contains parameters as defined in the configuration file 'congig.txt'
        
     Output:
        SeqDat_new:         dataframe; like input, but the original measurements are overwritten with separate samples from statistics
    '''
    import pandas as pd
    import numpy as np
    
    SeqDat_new = pd.DataFrame()
    ID_Col_Name = Name_Dict['ID_Col_Name']
    Sequence_column = Name_Dict['Sequence_column']
    Y_Col_Name = eval(Name_Dict['Y_Col_Name'])
    Stats_Std = eval(Name_Dict['Stats_Std'])
    for index, row in SeqDat.iterrows():
        myID = row[ID_Col_Name]
        mySeq = row[Sequence_column]
        Sample_number = min(row[eval(Name_Dict['Stats_Samples'])])
        Target_mean = [row[Y_Name] for Y_Name in Y_Col_Name]
        Target_std = [row[Std] for Std in Stats_Std]
        Sample_array = [stats2samples(mymean, mystd, mysamples) for mymean, mystd, mysamples in zip(Target_mean, Target_std, np.tile(Sample_number,(len(Y_Col_Name),1)))]
        for MySample in zip(*Sample_array):
            data = zip(list([ID_Col_Name]) + list([Sequence_column]) + Y_Col_Name, np.hstack([myID,mySeq,MySample]))
            SeqDat_new = SeqDat_new.append(dict(data), ignore_index=True)

    # correcting data types
    dtypes = dict(zip(Y_Col_Name, ['float64']*len(Y_Col_Name)))
    SeqDat_new = SeqDat_new.astype(dtypes)
    
    return SeqDat_new

def stats2samples(Target_mean, Target_std, Target_num=2):
    '''
    synthetic data generation based on information on mean, standard deviation, and sample size samples generated based on normal random number distribution 
    reference here: 
    https://stackoverflow.com/questions/51515423/generate-sample-data-with-an-exact-mean-and-standard-deviation
    
    Input:   
        Target_mean:    float; arithmetic mean of samples
        Target_std:     float; standard deviation of samples
        Target_num:     integer; number of samples to generate, default 2
    '''
    import numpy as np
    
    myx = np.random.normal(loc=Target_mean, scale=Target_std, size=Target_num)
    zeromean_data = (myx - np.mean(myx))
    zeromean_mean = np.mean(zeromean_data)
    zeromean_std = np.std(zeromean_data)
    scaled_data = zeromean_data * (Target_std/zeromean_std)
    scaled_mean = np.mean(scaled_data)
    scaled_std = np.std(scaled_data)
    final_data = scaled_data + Target_mean
#     final_mean = np.mean(final_data)
#     final_std = np.std(final_data)
    
    return final_data
    # 

###########################################################################
###########################################################################

def split_train_test(SeqDat, test_ratio=.1):
    '''
    Data split into training and test sets, with given ratio (default 10% test)
    
    Input:
        SeqDat:         dataframe; sequence in letter, label, and one-hot format, GC-content, expression strength and sequence IDs
        test_ratio: float [0,1), ratio of test data to total data
    
    Output:
        Arg#1:      dataframe; training data
        Arg#2:      dataframe; test data        

    '''
    import random 
    import numpy as np
    
    My_Observed,_ = SeqDat.shape
    Test_Idx = np.sort(np.array(random.sample(range(My_Observed), np.int(test_ratio* My_Observed))))
    Train_Idx = np.setdiff1d(range(My_Observed), Test_Idx)
    return SeqDat.iloc[Train_Idx], SeqDat.iloc[Test_Idx]

    
###########################################################################
###########################################################################
#
# Data analysis and synthetic library
#
###########################################################################
###########################################################################

# General function to identify conserved sequences

def Sequence_Conserved_Adjusted(SeqDat, Name_Dict, Entropy_cutoff=0):
    '''
    Adjusting sequence length by removal of positions with low variability
    
    Input:
        SeqDat:         dataframe; sequence in letter, label, and one-hot format, GC-content, expression strength and sequence IDs
        Name_Dict:      dictionary; contains parameters as defined in the configuration file 'congig.txt'
        Entropy_cutoff: float; threshold for variability report in entropy values, default=0
        
    Output:
        SeqDat:          dataframe; like input with sequences shortened by non-informative positions
    '''
    import numpy as np
    
    # Position_Conserved represents position for which only a single nucleotide was measured, and thus no diversity/information exists
    Position_Conserved, PSEntropy = Conserved_Sequence_Exclusion(np.array(SeqDat['Sequence_label-encrypted'].tolist()), Entropy_cutoff)
    # if there are positions without information in Position_Conserved they are deleted
    if bool(Position_Conserved.size):
        SeqDat = SeqDat.assign(ColName=SeqDat['Sequence_label-encrypted'])
        SeqDat.rename(columns={'ColName':'Sequence_label-encrypted_full'}, inplace=True);
        SeqDat['Sequence_label-encrypted'] = list(np.delete(np.array(list(SeqDat['Sequence_label-encrypted'])),Position_Conserved, axis=1))
    SeqDat['OneHot'] = list_onehot(SeqDat['Sequence_label-encrypted'])
    
    return SeqDat, Position_Conserved, PSEntropy

def Conserved_Sequence_Exclusion(SeqLab, Entropy_cutoff=0):
    '''
    Returns the sequence positions that have a lower or equal diversity than given as entropy threshold.
    
    Input:
          SeqLab:              np-array; columns represent sequence position, rows represent samples
          Entropy_cutoff:      float; threshold for variability report in entropy values, default=0
          
    Output:
          Position_Conserved:  np-array; positions of low variability within the sequence length
     '''
#    import pandas as pd
    import numpy as np
    
    PSEntropy = Entropy_on_Position(SeqLab)
    Position_Conserved = np.arange(len(PSEntropy))[PSEntropy <= Entropy_cutoff]
    
    return Position_Conserved, PSEntropy

###########################################################################
###########################################################################

def Sequence_Dist_DiffSum(SeqObj):
    '''
    Returns the genetic sequence distance all sequences in the data list. The distance is determined from the sum of difference in bases divided by total base number, i.e. max difference is 1, identical sequence is 0.
    
    Input:
         SeqObj:             list; From original Data, the sequence in conventional letter format
         
    Output:
         PromDist_SymMatrix: np-array; genetic distances
    '''
    import numpy as np

    Num_Samp = len(SeqObj)
    PromDist = list()
    for idx1 in range(Num_Samp):
        for idx2 in range(idx1, Num_Samp):
            PromDist.append(np.sum([int(seq1!=seq2) for seq1,seq2 in zip(SeqObj[idx1],SeqObj[idx2])], dtype='float')/len(SeqObj[idx1]))
    
    Entry_Size = np.insert(np.cumsum(np.arange(Num_Samp,0,-1)),0,0)
    PromDist_SymMatrix = np.zeros((Num_Samp,Num_Samp))
    for index in range(0,Num_Samp):
        PromDist_SymMatrix[index,index:] = PromDist[Entry_Size[index]:Entry_Size[index+1]]
        PromDist_SymMatrix[index:,index] = PromDist[Entry_Size[index]:Entry_Size[index+1]]

    return PromDist_SymMatrix

###########################################################################
###########################################################################

def Sequence_Ref_DiffSum(SeqObj):
    '''
    Returns the genetic sequence distance relative to the first sequence in the data-frame to all following sequences in the data list. The distance is determined from the sum of difference in bases divided by total base number, i.e. max difference is 1, identical sequence is 0.
    
    Input:
           SeqObj: list; From original Data, the sequence in conventional letter format
           
    Output:
           PromDist: np-array; genetic distances 
    '''
    import numpy as np

    RefSeq = SeqObj[0]
    Num_Samp = len(SeqObj)
    PromDist = list()
    for idx1 in range(1, Num_Samp):
        PromDist.append(np.sum([int(seq1!=seq2) for seq1,seq2 in zip(RefSeq, SeqObj[idx1])], dtype='float')/len(SeqObj[idx1]))
    
    return np.array(PromDist)

###########################################################################
###########################################################################

def Find_Near_Seq(SeqTest, SeqRef):
    '''
    Identifies the closes Reference sequence to a test sequence
    
    Input:      SeqTest:      string; single sequence
                SeqRef:       list; array
    
    Output:     SeqTarget:    string; sequence from reference with closest distance to test
                SeqClost_Idx: integer; id of the closes sequence to test
    '''
    import numpy as np
    
    # irrelevant sequence positions can be marked by 'N', here we find them and only use true bases for finding the target sequence
    PosTrue = [i for i,x in enumerate(SeqTest) if x != 'N']
    Num_Samp = len(SeqRef)
    PromDist = list()
    for idx1 in range(Num_Samp):
        # only selecting relevant bases
        seqtest = [SeqTest[i] for i in PosTrue]
        seqref = [SeqRef[idx1][i] for i in PosTrue]
        PromDist.append(np.sum([int(seqtest!=seqref) for seqtest,seqref in zip(seqtest, seqref)], dtype='float')/len(PosTrue))
    
    PromDist = np.array(PromDist)
    SeqClose_Idx = np.argsort(PromDist)[0]
    
    SeqTarget = np.unique(SeqRef[SeqClose_Idx])
    
    return SeqTarget, SeqClose_Idx

###########################################################################
###########################################################################

def Extract_MultLibSeq(SeqDat, Target, Seq_Numb, Y_Col_Name):
    '''
    Finding the closest measured samples in experimental library to the target expression 
    
    Input:
        SeqDat:      dataframe; sequence in letter, label, and one-hot format, GC-content, expression strength and sequence IDs
        Target:      np-array; target expression strength
        Y_Col_Name:  string; name of expression strength column
      
    Output:
        SeqObj:      list; sequences close to target expression
        Target_lst:  list; indices of close sequences to target expression
    
    '''
    import numpy as np
    
    Target_lst = []
    SeqObj = []
    Expr_Prom_Num = Target.shape[1]
    
    for Expr_idx in range(Expr_Prom_Num):
        # finding positions of unique sequences
        u, index = np.unique(SeqDat['Sequence_letter-encrypted'].str.upper(), return_inverse=True)
        Base_Expr = SeqDat[Y_Col_Name].values
        # calculation of the distance of each promoter strength to target strength
        TR_Dist = np.abs(1-Base_Expr/Target[:,Expr_idx])
        # sorting the closest expression strength to the beginning
        Target_Idx = np.argsort(np.sum(TR_Dist, axis=1))
        # replicates should all have similar expression values, but we want to know different sequences that are close to the target sequence
        # We insert the ordered expression distance into the categories from the unique sequences
        # the output gives the position of the closest expression in the categorized sequences
        _, i2 = np.unique(index[Target_Idx], return_index=True)
        Seq_lst = u[np.argsort(i2)[:Seq_Numb]]

        # The extraction of the position of the closest expression in the original data is complex because of replicates.
        # In the following we select the closest indices for unique sequences
        Expr_ord = index[Target_Idx]
        Expr_uni = []
        [Expr_uni.append(x) for x in Expr_ord if x not in Expr_uni]  
        Idx_lst = np.array([list(Expr_ord).index(Indx) for Indx in Expr_uni])
        Target_unique = Target_Idx[Idx_lst]

        Target_lst.append(Target_unique[:Seq_Numb])
        SeqObj.append(Seq_lst)
    
    return SeqObj, Target_lst

###########################################################################
###########################################################################

# Entropy calculation
def Entropy_on_Position(PSArray):
    '''
    Analysis of position specific entropy. 
    
    Input: 
        PSArray:   np array; columns represent sequence position, rows represent samples
        
    Output:
        PSEntropy: np vector; entropy of each sequence position
    '''
    import numpy as np
    from scipy.stats import entropy
    
    PSEntropy = list()
    for col in PSArray.T:
        value, counts = np.unique(col, return_counts=True)
        PSEntropy.append(entropy(counts, base=2))
        
    return np.array(PSEntropy)

###########################################################################
###########################################################################
#
# Regression functions
#
###########################################################################
###########################################################################

def MyRF(SeqOH, Validation_cutoff=.1, Num=100, Y_Col_Name='promoter activity', Response_Value=1, AddFeat=None):
    '''
    This function trains a random forest regressor and performs gradient search for optimal parameters with group shuffle shift.
    
    Input:
        SeqOH:             dataframe; sequence in letter, label, and one-hot format, GC-content, expression strength and sequence IDs
        Validation_cutoff: float; ratio of cross-validation train and test sets
        Num:               integer; number of validation splits performed
        Y_Col_Name:        string; name of expression strength column
        Response_Value:    integer; decision how the original expression values where manipulated, 0: standardized, 1: original, >1: categorization
        AddFeat:           string; name of additional feature, typically GC-content
        
    Output:
        grid_forest:       function; regressor 
    '''
    from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
    from sklearn.model_selection import GroupShuffleSplit, GridSearchCV
    import numpy as np
       
    Sequence_Samples, Sequence_Positions, Sequence_Bases = np.array(SeqOH['OneHot'].values.tolist()).shape
    X = np.array(SeqOH['OneHot'].values.tolist()).reshape(Sequence_Samples,Sequence_Positions*Sequence_Bases)
    # adding rows to x for additional features
    if AddFeat != None:
        X = np.append(X,np.array(SeqOH[AddFeat]), axis=1)
    Y = SeqOH[Y_Col_Name].values
    groups = SeqOH['Sequence_letter-encrypted']
    Number_Estimators = np.arange(20,40,2)
    Max_Features = np.arange(10,30,2)
    min_samples_split = np.arange(3,4,1)
    max_depth = np.arange(3,5)
    min_samples_leaf = np.array([2])
    param_grid = [{'bootstrap':[False], 'n_estimators': Number_Estimators, 'max_features': Max_Features, 'min_samples_split': min_samples_split, 'max_depth':max_depth, 'min_samples_leaf':min_samples_leaf}]
    # Group shuffle split removes groups with identical sequences from the development set
    # This is more realistic for parameter estimation
    cv = GroupShuffleSplit(n_splits=Num, test_size=Validation_cutoff, random_state=42)

    if Response_Value > 1:
        forest_grid = RandomForestClassifier()
    else:
        forest_grid = RandomForestRegressor()
    grid_forest = GridSearchCV(forest_grid, param_grid, cv=cv, n_jobs=-1)
    grid_forest.fit(X, Y, groups)
           
    return grid_forest

###########################################################################
###########################################################################

def MyGB(SeqOH, Validation_cutoff=.1, Num=100, Y_Col_Name='promoter activity', Response_Value=1, AddFeat=None):
    '''
    This function trains a gradient boosting regressor and performs gradient search for optimal parameters with group shuffle shift.
    
    Input:
        SeqOH:             dataframe; sequence in letter, label, and one-hot format, GC-content, expression strength and sequence IDs
        Validation_cutoff: float; ratio of cross-validation train and test sets
        Num:               integer; number of validation splits performed
        Y_Col_Name:        string; name of expression strength column
        Response_Value:    integer; decision how the original expression values where manipulated, 0: standardized, 1: original, >1: categorization
        AddFeat:           string; name of additional feature, typically GC-content
        
    Output:
        grid_GBR:       function; regressor 
    '''
    from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier
    from sklearn.model_selection import GroupShuffleSplit, GridSearchCV
    import numpy as np

    Sequence_Samples, Sequence_Positions, Sequence_Bases = np.array(SeqOH['OneHot'].values.tolist()).shape
    X = np.array(SeqOH['OneHot'].values.tolist()).reshape(Sequence_Samples,Sequence_Positions*Sequence_Bases)
    # adding rows to x for additional features
    if AddFeat != None:
        X = np.append(X,np.array(SeqOH[AddFeat]), axis=1)
    Y = SeqOH[Y_Col_Name].values

    groups = SeqOH['Sequence_letter-encrypted']
    Number_Estimators = np.arange(20,40,3)
    Max_Features = np.arange(10,30,3)
    min_samples_split = np.arange(2,4,1)
    learning_rate = np.logspace(-3,2,5)
    param_grid = [{'n_estimators': Number_Estimators, 'max_features': Max_Features, 'min_samples_split': min_samples_split, 'learning_rate': learning_rate}]
    # Group shuffle split removes groups with identical sequences from the development set
    # This is more realistic for parameter estimation
    cv = GroupShuffleSplit(n_splits=Num, test_size=Validation_cutoff, random_state=42)

    if Response_Value > 1:
        GB = GradientBoostingClassifier()
    else:
        GB = GradientBoostingRegressor()
    grid_GB = GridSearchCV(GB, param_grid, cv=cv, n_jobs=-1)
    grid_GB.fit(X, Y, groups)   
    
    return grid_GB

###########################################################################
###########################################################################

def MySV(SeqOH, Validation_cutoff=.1, Num=100, Y_Col_Name='promoter activity', Response_Value=1, AddFeat=None):
    '''
    This function trains a support vector regressor and performs gradient search for optimal parameters with group shuffle shift.
    
    Input:
        SeqOH:             dataframe; sequence in letter, label, and one-hot format, GC-content, expression strength and sequence IDs
        Validation_cutoff: float; ratio of cross-validation train and test sets
        Num:               integer; number of validation splits performed
        Y_Col_Name:        string; name of expression strength column
        Response_Value:    integer; decision how the original expression values where manipulated, 0: standardized, 1: original, >1: categorization
        AddFeat:           string; name of additional feature, typically GC-content
        
    Output:
        grid_SVR:       function; regressor 
    '''
    from sklearn import svm
    from sklearn.model_selection import GroupShuffleSplit, GridSearchCV
    import numpy as np

    Sequence_Samples, Sequence_Positions, Sequence_Bases = np.array(SeqOH['OneHot'].values.tolist()).shape
    X = np.array(SeqOH['OneHot'].values.tolist()).reshape(Sequence_Samples,Sequence_Positions*Sequence_Bases)
    # adding rows to x for additional features
    if AddFeat != None:
        X = np.append(X,np.array(SeqOH[AddFeat]), axis=1)
    Y = SeqOH[Y_Col_Name].values

    groups = SeqOH['Sequence_letter-encrypted']
    C_values = np.logspace(-3,3,20)
    gamma_values = np.logspace(-3,1,20)
    param_grid = [{'C': C_values, 'gamma': gamma_values, 'kernel': ['rbf']}]
    # Group shuffle split removes groups with identical sequences from the development set
    # This is more realistic for parameter estimation
    cv = GroupShuffleSplit(n_splits=Num, test_size=Validation_cutoff, random_state=42)

    if Response_Value > 1:
        SV = svm.SVC()
    else:
        SV = svm.SVR()
    grid_SV = GridSearchCV(SV, param_grid, cv=cv, n_jobs=-1)
    grid_SV.fit(X, Y, groups)     
    
    return grid_SV

###########################################################################
###########################################################################

def SequenceRandomizer_Parallel(RefSeq, Base_SequencePosition, n=1000):
    '''
    This function generates random sequence combinations. It takes the reference sequence and changes nucleotides at positions that have been experimentally tested. Only as much nucleotides are changed to remain within a given sequence distance.
    '''
    import numpy as np
    import multiprocessing
    from joblib import Parallel, delayed
#     from ExpressionExpert_Functions import SequenceRandomizer_Single
    
    num_cores = multiprocessing.cpu_count()
    use_core = min([num_cores, n])
  
    Result = Parallel(n_jobs=use_core)(delayed(SequenceRandomizer_Single)(RefSeq, Base_SequencePosition) for idx in range(n))

    return Result # Sequence_multiple

###########################################################################
###########################################################################

def SequenceRandomizer_Single(RefSeq, Base_SequencePosition):
    
    import numpy as np
    import random
    from ExpressionExpert_Functions import list_integer
    
    Alphabet = ['A','C','G','T']
    
    #   Maximum number of nucleotides that can be changed simultaneously based on the sequence distance cut-off
    Nucleotide_Replace_Numb = len(Base_SequencePosition)
    # Generating the positions with nucleotide randomization
    MySynChange = np.array(random.sample(list(Base_SequencePosition.index), Nucleotide_Replace_Numb))
    # the following dataframe has the information of experimentally tested nucleotides as boolean table
    mytmp = Base_SequencePosition.loc[MySynChange]
    # The boolean table of tested nucleotides is converted into an array containing the explicit nucleotide letters
    myArr = np.tile(Alphabet, (Nucleotide_Replace_Numb,1))
    # following, non-tested nucleotides are deleted
    Pos_Del, Nucl_Del = np.where(mytmp.values == 0)
    if not Pos_Del.any(): # != '':
        # deleting non-tested nucleotides
        myArr[tuple([Pos_Del,Nucl_Del])] = 'X'
    # Generating a reference sequence to work with
    TstSeq = list(RefSeq)
    # Changing indices from nucleotide oriented to array oriented
    ArSynChange = MySynChange + len(RefSeq)
    # setting the reference nucleotides at the positions to 'X', if they should not be chosen for new random sequences. This maximizes the distance of all random sequences, dont use it by default
    # RefSeq_IntLab = np.array(list_integer(RefSeq))
    # myArr[[range(Nucleotide_Replace_Numb)], [np.reshape(RefSeq_IntLab[ArSynChange], -1)]] = 'X'

    # converting the nucleotide array to a list, so we can delete non-tested nucleotides
    Position_list = myArr.tolist()
    Seq_Base = list()
    for Position in Position_list:
        Seq_Base.append(list(set(Position).difference(set('X'))))
    # print(Seq_Base)
    # randomly choosing a possible nucleotide over the total number of exchange positions
    Replace_Bases = [PosIdx[np.random.randint(len(PosIdx))] for PosIdx in Seq_Base]
    #   Replacing the bases in the reference sequence
    for MutIdx, MutNucl in zip(ArSynChange, Replace_Bases):
        TstSeq[MutIdx] = MutNucl    
    Sequence_Single = ''.join(TstSeq)
    
    return Sequence_Single

###########################################################################
###########################################################################
# visualization
###########################################################################
###########################################################################

def ExpressionStrength_HeatMap(SeqList, Y_Col_Name = 'promoter activity'):
    '''
    Calculating the base and position specific average expression strength.
    Input:
        SeqList_df: dataframe, contains OneHot encoded sequence labelled 'OneHot' and 
                            expression strength parameters labelled as in variable 'Y_Col_Name'
    Output:
        Expression_HeatMap_df: dataframe, contains the average expression strength for each base on each position
    '''
    import numpy as np
    import pandas as pd
    
    # extracting the One-Hot encoding for all sequences
    my_rows = SeqList['OneHot'].shape
    Seq_OneHot_ar = np.array(list(SeqList['OneHot'])).reshape(my_rows[0],-1)
    # getting the number of positions considered
    myBasePos = Seq_OneHot_ar.shape[1]
    Exp_mult_ar = np.tile(np.array(SeqList[Y_Col_Name]),[myBasePos,1]).T
    Expr_OneHot_mean_ar = np.mean(np.multiply(Seq_OneHot_ar, Exp_mult_ar), axis=0)
    Expr_OneHot_mean_ar[Expr_OneHot_mean_ar==0.] = np.nan
    Expression_HeatMap_df = pd.DataFrame(Expr_OneHot_mean_ar.reshape(-1,4), columns=['A','C','G','T'])
    return Expression_HeatMap_df

###########################################################################
###########################################################################

def ExpressionVariation_HeatMap(SeqList, Y_Col_Name = 'promoter activity'):
    '''
    Calculating the base and position specific variation of expression strength.
    Input:
        SeqList_df: dataframe, contains OneHot encoded sequence labelled 'OneHot' and 
                            expression strength parameters labelled as in variable 'Y_Col_Name'
    Output:
        Variation_HeatMap_df: dataframe, contains the variation of expression strength for each base on each position
    '''
    import numpy as np
    import pandas as pd
    
    # extracting the One-Hot encoding for all sequences
    my_rows = SeqList['OneHot'].shape
    Seq_OneHot_ar = np.array(list(SeqList['OneHot'])).reshape(my_rows[0],-1)
    # getting the number of positions considered
    myBasePos = Seq_OneHot_ar.shape[1]
    Exp_mult_ar = np.tile(np.array(SeqList[Y_Col_Name]),[myBasePos,1]).T
    Expr_OneHot_var_ar = np.var(np.multiply(Seq_OneHot_ar, Exp_mult_ar), axis=0)
    Expr_OneHot_var_ar[Expr_OneHot_var_ar==0.] = np.nan
    Expression_HeatMap_df = pd.DataFrame(Expr_OneHot_var_ar.reshape(-1,4), columns=['A','C','G','T'])
    return Expression_HeatMap_df

###########################################################################
###########################################################################

def ExpressionMean_ttest(a_np, b_np, OH_List):
    '''
    Calculating whether the mean expression values are significantly different.
    Input:
        a_np:       np vector, contains the measurements of cross host/reporter
        b_np:       np vector, contains the measurements of cross host/reporter
        OH_List:    list, contains OneHot encoded sequence in list format generate by e.g. list(SeqDat['OneHot'])
    Output:
        Variation_HeatMap_df: dataframe, contains the variation of expression strength for each base on each position
    '''
    import numpy as np
    import pandas as pd
    from scipy.stats import ttest_ind
    
    # creating the filter
    Seq_OneHot_ar = np.array(OH_List).reshape(len(OH_List),-1)
    myBasePos = Seq_OneHot_ar.shape[1]
    # transforming the input to the filter dimensions
    Expr_mult_a = np.tile(a_np, [myBasePos,1]).T
    Expr_mult_b = np.tile(b_np, [myBasePos,1]).T
    # multiplying the input with the filter
    a_matr = np.multiply(Seq_OneHot_ar, Expr_mult_a)
    b_matr = np.multiply(Seq_OneHot_ar, Expr_mult_b)

    # performing the students t-test
    stat, p = ttest_ind(a_matr, b_matr, axis=0)
    # identifying significant differences in the mean values
    MeanSignDiff = np.multiply(p<=.05,1).reshape(-1,4)
    Expression_ttest = pd.DataFrame(MeanSignDiff, columns=['A','C','G','T'])
    return Expression_ttest

###########################################################################
###########################################################################

def df_HeatMaps(Data_df, Z_Label, Plot_Save=False, Plot_File='dummy', cbar_lab=None, FigFontSize=14):
    '''
    Function for heat map generation
    Input:
            Data_df:     dataframe, rows represent sequence position, columns represent bases with their names as labels
            Z_Label:     string, label for the colorbar
            Plot_Save:   boolean, decision whether to save plot
            Plot_File:   string, name for the figure file
            cbar_lab:    sting vector, names for the color bar
            FigFontSize: integer, font size for axis labels
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig, ax = plt.subplots()
    im = ax.pcolor(Data_df.T, cmap='rainbow')
        
    # identifying dataframe dimensions
    my_rows, my_cols = Data_df.shape
    
    #label names
    my_Seq_range = np.arange(0,my_rows,5)
    row_labels = np.arange(-my_rows,0,5)
    col_labels = Data_df.T.index
    #move ticks and labels to the center
    ax.set_xticks(ticks=my_Seq_range+0.5, minor=False)
    ax.set_yticks(np.arange(Data_df.T.shape[0])+0.5, minor=False)
    #insert labels
    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(col_labels, minor=False)
    #rotate label if too long
#     plt.xticks(rotation=90)

    plt.xlabel('Position', fontsize=FigFontSize)
    plt.ylabel('Base', fontsize=FigFontSize)
    # plt.title('asdf')
    if cbar_lab is not None:
        tick_num = len(cbar_lab)
        myDat = Data_df.values.reshape(-1)
        myticks = np.histogram(myDat[~np.isnan(myDat)], bins=tick_num-1)[1]
        cbar = fig.colorbar(im, ticks=myticks, label=Z_Label)
        cbar.ax.set_yticklabels(cbar_lab)
    else:
        fig.colorbar(im, label=Z_Label) # , fontsize=FigFontSize
        
    if Plot_Save:
        Fig_Type = Plot_File.split('.')[1]
        plt.savefig(Plot_File, format=Fig_Type)
    plt.show()

###########################################################################
###########################################################################
       
def generate_distances(SeqDat, Name_Dict, Hist_Type):
    '''
    Function to generate sequence distances.
    
    Input:
        SeqDat:        dataframe, contains sequences
        Name_Dict:     dictionary, potentially contains reference sequence
        Hist_Type:     boolean, decision to determine sequence distance relative to all sequences (0) or to a reference sequence (1)
    
    Output:
        mydist:        matrix, contains pairwise sequence distances, either to reference or to all sequences
    '''
    import numpy as np
    
    # if there are too many samples pair-wise distance over all samples is too costly computationally
    # The distance based on reference sequence is similar to the pairwise distance
    if len(SeqDat) > 1000:
        Hist_Type = 1
    
    # distances are determined relative to reference
    if Hist_Type == 1:
        if Name_Dict['RefSeq'] is not '':
            RefSeq = Name_Dict['RefSeq']
            print('use reference sequence')
        else:    
            # using the one-hot encoding the most common nucleotide on each position is calculated.
            Alphabet = ['A','C','G','T']
            Pos_Nucl_Sum = np.sum(np.dstack(SeqDat['Sequence'].values), axis=2)
            RefSeq_list = list([Alphabet[Pos_idx] for Pos_idx in np.argmax(Pos_Nucl_Sum, axis=1)])
            RefSeq = ''.join(RefSeq_list)
    
        print('Distance reference sequence:', RefSeq)
        SeqDat_wRef = list(SeqDat['Sequence_letter-encrypted'])
        SeqDat_wRef.insert(0, RefSeq)

        mydist = Sequence_Ref_DiffSum(SeqDat_wRef)
    # distances are calculated among all sequences
    else: 
        print('Complete sequence distances.')
        dist_tmp = Sequence_Dist_DiffSum(SeqDat['Sequence'].values)
        mydist = dist_tmp[np.triu_indices(dist_tmp.shape[1], k=1)]
    
    return mydist

###########################################################################
###########################################################################
    
def my_CrossValScore(X, Y, groups, cv, ML_fun, metric):
    '''
    Function to generate statistics according to cross validation scoring
    '''
    import numpy as np
    
    mytrain = list()
    mytest = list()
    tst_num = cv.n_splits
    for idx in range(tst_num):
        cv.n_splits = 1
        mycv = list(cv.split(X, Y, groups))
        onetrain = mycv[0][0]
        onetest = mycv[0][1]
        ML_fun.fit(X[onetrain], Y[onetrain])
        mytrain.append({'train':metric(ML_fun, X[onetest], Y[onetest])})
#         myscores.append(metric(ML_fun, X[onetest], Y[onetest]))
        
    myscores = dict({'TrainR2':mytrain})
    return myscores
    
###########################################################################
###########################################################################
    
def DetectOutlier(myData):
    '''
    Function to identify outlier in one dimensional input. Outlier are detected as being more distant than three times the standard deviation from the mean value.
    '''
    import numpy as np
    
    myMean = np.mean(myData)
    myStd = np.std(myData)

    # Outliers defined as exceeding 3x std
    OutThresh = 3
    LowerLimit = myMean - myStd*OutThresh
    UpperLimit = myMean + myStd*OutThresh

    myOutlier = [indx for indx, val in enumerate(myData) if val > UpperLimit or val < LowerLimit]

    return myOutlier
    
###########################################################################
###########################################################################
    
def Predict_SequenceActivity(Sequence, Name_Dict):
    '''
    Function to predict expression activity for single input sequence. The embedding of the input sequence in the exploratory space is tested with the sequence distance but not position diversity.
    
    Input:
        Sequence:    string, sequence
        Name_Dict:   dictionary, config information
        
    Output:
        Activity:    float, expression activity for each promoter library
    '''
    from sklearn.preprocessing import OneHotEncoder
    import numpy as np
    import pandas as pd
    import os
    import joblib
    import pickle

    Measure_Numb = int(Name_Dict['Library_Expression'])
    Data_File = Name_Dict['Data_File']
    # extract the filename for naming of newly generated files
    File_Base = Data_File.split('.')[0]
    # the generated files will be stored in a subfolder with custom name
    Data_Folder = 'data-{}'.format(File_Base)
    Y_Col_Name = eval(Name_Dict['Y_Col_Name'])
    ML_Date = Name_Dict['ML_Date']
    ML_Type = Name_Dict['ML_Regressor']
    encoder = OneHotEncoder()
    SeqOH = encoder.fit_transform(np.array(list_integer(Sequence)))
    
    # calculating GC content
    GC_content = (Sequence.count('G')+Sequence.count('C'))/len(Sequence)
    
    Activity = np.empty(Measure_Numb)
    # testing whether sequence is within exploratory space
    SeqDF = pd.DataFrame({'Sequence_letter-encrypted': [Sequence], 'Sequence': [SeqOH.toarray()]})
    mydist = generate_distances(SeqDF, Name_Dict, Hist_Type=1)
    if mydist < eval(Name_Dict['Sequence_Distance_cutoff']):
        ML_Best = dict()
        Expr_Scaler = dict()
        for Meas_Idx in range(Measure_Numb): 
        #     print('loading Regressor #{}'.format(Meas_Idx))
            Regressor_File = os.path.join(Data_Folder, '{}_{}_{}_{}-Regressor.pkl'.format(ML_Date, File_Base, Y_Col_Name[Meas_Idx].replace(' ','-'), ML_Type))
            Parameter_File = os.path.join(Data_Folder, '{}_{}_{}_{}-Params.pkl'.format(ML_Date, File_Base, Y_Col_Name[Meas_Idx].replace(' ','-'), ML_Type))
            ML_DictName = '{}_Regressor'.format(Y_Col_Name[Meas_Idx])
            ML_Best[ML_DictName] = joblib.load(Regressor_File)
            # I assume the parameters have been generated in the same run as the regressor itself and is located in the same directory following the default naming scheme
            Data_Prep_Params = pickle.load(open(Parameter_File,'rb'))
            # extracting the positions that were removed because of insufficient information content
            Positions_removed = Data_Prep_Params['Positions_removed']
            # removed positions should be identical to the reference
            if Name_Dict['RefSeq'] is not '':
                Sequence_relevant = ''.join((char for idx, char in enumerate(Sequence) if idx in Positions_removed))
                Reference_relevant = ''.join((char for idx, char in enumerate(Name_Dict['RefSeq']) if idx in Positions_removed))
                if Sequence_relevant != Reference_relevant:
                    print('Mutations out of exploratory space')
                    return
            # if the data was standardized we load the corresponding function
            if eval(Name_Dict['Data_Standard']):
                # extracting standard scaler from existing random forest regressor
                # The standard scaler default name is the name of the expression measurement column with suffix: '_Scaler'
                Scaler_DictName = '{}_Scaler'.format(Y_Col_Name[Meas_Idx])
                Expr_Scaler[Scaler_DictName] = Data_Prep_Params[Scaler_DictName]
            X_tmp = np.hstack(np.delete(SeqOH.toarray(), Positions_removed, axis=0))
            # adding overall GC content
            X = np.append(X_tmp, GC_content).reshape(1,-1)
            Activity[Meas_Idx] = ML_Best[ML_DictName].predict(X)
            # correcting the prediction for standardized data
            if eval(Name_Dict['Data_Standard']):
                Activity[Meas_Idx] = Expr_Scaler[Scaler_DictName].inverse_transform(Activity[Meas_Idx].reshape(1,-1))
    else:
        print('Input sequence is too far from the reference sequence and cannot be reliably predicted.')
        return
    
    return np.round(Activity,3)