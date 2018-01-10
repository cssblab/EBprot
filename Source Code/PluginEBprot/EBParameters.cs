using BaseLibS.Num;
using BaseLibS.Param;
using PerseusApi.Matrix;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PluginEBprot
{
    public static class EBParameters
    {
        public const int EBPROTMODE = 0;
        public const int GRPMODE = 1;
        public const string DIR = "Working Directory";
        public const string PEPID = "Peptide Identifier Column";
        public const string PROTID = "Protein Identifier Column";
        public const string LABELNUM = "# of Labels";
        public const string GROUPNUM = "# of Groups";
        public const string MINPEP = "Min # of Peptides";
        public const string DATAFORM = "Data Input Form";
        public const string EXPDESIGN = "Experimental Design";
        public const string OUTFILT = "Outlier Filtering";
        public const string MINK = "Min K";
        public const string BFDRCUT = "Bayesian False Discovery Rate Cutoff";
        public const string LEFTB = "Lower Bound";
        public const string RIGHTB = "Upper Bound";
        public const string CONTRAST = "Contrast";
        public const string MINSAMP = "Min # of Samples";
        public static FolderParam GetWorkingDir()
        {
            return new FolderParam(DIR)
            {
                Help = "Specifies the directory where internal input files, output matrices and plots produced by EBprot will be saved.",
            };
        }

        public static SingleChoiceParam GetPeptideIDColumn(IMatrixData mdata)
        {
            return new SingleChoiceParam(PEPID, 0)
            {
                Values = mdata.StringColumnNames,
                Help = "The peptide identifiers in EBprot analysis."
            };
        }

        public static SingleChoiceParam GetProteinIDColumn(IMatrixData mdata)
        {
            return new SingleChoiceParam(PROTID, 1)
            {
                Values = mdata.StringColumnNames,
                Help = "The protein identifiers in EBprot analysis."
            };
        }

        //Version 1 
        //public static MultiChoiceParam GetData(IMatrixData mdata, string name, string help)
        //{
        //    string[] valuesList = ArrayUtils.Concat(mdata.ColumnNames, mdata.NumericColumnNames);
        //    int[] def = valuesList.Length > 0 ? Enumerable.Range(0, valuesList.Length).ToArray() : new int[0];
        //    return new MultiChoiceParam(name, def)
        //    {
        //        Values = valuesList,
        //        Repeats = false,
        //        Help = help
        //    };
        //}

        public static MultiChoiceParam GetData(IMatrixData mdata, string name, string help, int[] def)
        {
            string[] valuesList = ArrayUtils.Concat(mdata.ColumnNames, mdata.NumericColumnNames);
            return new MultiChoiceParam(name, def)
            {
                Values = valuesList,
                Repeats = false,
                Help = help
            };
        }

        //Version 1 
        //public static Parameter[] SelectData(IMatrixData mdata, int mode)
        //{
        //    string name = null;
        //    string help = null;
        //    //fill
        //    if (mode == EBProtMode)
        //    {
        //        name = "Intensity Data";
        //        help = "";
        //    }
        //    else if (mode == GRPMode)
        //    {
        //        name = "Intensity Data";
        //        help = "";
        //    }

        //    return new Parameter[]
        //    {
        //        GetPeptideIDColumn(mdata),
        //        GetProteinIDColumn(mdata),
        //        GetData(mdata, name, help)
        //    };
        //}

        public static string LabelString(int i, string str)
        {
            return "Label " + i.ToString() + " " + str;
        }

        public static SingleChoiceWithSubParams GetLabels(IMatrixData mdata)
        {
            string[] choiceList = new string[10];
            //be careful about the difference between Parameter and Parameters!
            Parameters[] subParamList = new Parameters[10];
            int parNameWidth = 100;
            int totWidth = 600;

            for (int i = 0; i < 10; i++)
            {
                choiceList[i] = (i + 1).ToString();
                Parameter[] subsubParam = new Parameter[2*(i + 1)];
                //how many data for each label when equally split
                int split = ArrayUtils.Concat(mdata.ColumnNames, mdata.NumericColumnNames).Length / (i + 1);
                for (int j = 0; j <= i; j++)
                {
                    subsubParam[2 * j] = new StringParam(LabelString(j + 1, "Name"));
                    int[] def = split != 0 ? Enumerable.Range(j*split, split).ToArray() : new int[0];
                    //fill
                    subsubParam[2 * j + 1] = GetData(mdata, LabelString(j + 1, "Data"), "", def);
                }
                subParamList[i] = new Parameters(subsubParam);

            }

            //Parameters param1 = new Parameters(new Parameter[] { new StringParam("Label 1 Name"), GetData(mdata, "Label 1 Data", "") });
            //Parameters param2 = new Parameters(new Parameter[] { new StringParam("Label 1 Name"), GetData(mdata, "Label 1 Data", ""), new StringParam("Label 2 Name"), GetData(mdata, "Label 2 Data", "") });
            //string[] choiceList = new[] { "1", "2" };
            //Parameters[] subParamList = new[] { param1, param2 };
            return new SingleChoiceWithSubParams(LABELNUM, 0)
            {
                Values = choiceList,
                SubParams = subParamList,
                ParamNameWidth = parNameWidth,
                TotalWidth = totWidth,
                //fill
                Help = ""
            };
        }

        public static IntParam GetMinPep()
        {
            return new IntParam(MINPEP, 1)
            {
                Help = "Minimum number of peptides per protein required in the analysis."
            };
        }

        public static Parameter[] SelectEBData(IMatrixData mdata)
        {
            //fill
            //string name = "Intensity Data";
            //string help = "";

            return new Parameter[]
            {
                GetMinPep(),
                GetPeptideIDColumn(mdata),
                GetProteinIDColumn(mdata),
                //GetData(mdata, name, help)
                GetLabels(mdata),
            };
        }

        public static SingleChoiceWithSubParams GetDataForm()
        {
            int parNameWidth = 100;
            int totWidth = 600;

            Parameters emptySubParams = new Parameters();
            Parameters customSubParamsBase = new Parameters(new Parameter[] { new DoubleParam("Base", 2.0) });
            string[] baseList = new[] { "Raw", "log_2", "ln", "log_10", "log_custom" };
            Parameters[] subParamListBase = new[] { emptySubParams, emptySubParams, emptySubParams, emptySubParams, customSubParamsBase };

            return new SingleChoiceWithSubParams(DATAFORM, 0)
            {
                Values = baseList,
                SubParams = subParamListBase,
                ParamNameWidth = parNameWidth,
                TotalWidth = totWidth,
                Help = "Specification for the data input form.\nRaw: unprocessed, untransformed data.\nln: log_e transformed data.\nlog_2: log_2 transformed data.\nlog_10:  log_10 transformed data.\nlog_custom: log_X transformed data, where X is a specified positive real value\nIf not log2 then will be converted to that."
            };
        }

        public static SingleChoiceParam GetExpDesign()
        {
            return new SingleChoiceParam(EXPDESIGN, 0)
            {
                Values = new [] { "Independent", "Replicate", "Timecourse" },
                //fill
                Help = ""
            };
        }

        public static Parameter[] GetAboutEBData()
        {
            return new Parameter[]
            {
                GetExpDesign(),
                GetDataForm()
            };
        }

        public static BoolWithSubParams GetOutlierRM() => new BoolWithSubParams(OUTFILT, true)
        {
            Default = true,
            Help = "Whether or not to perform the outlier-filtering procedure. If true, peptides with aberrant ratios (falling outside of the reference stan- dard deviation) from each protein that has at least 3 ratios will be removed. Otherwise, the data will be analysed as it is.",
            SubParamsTrue = new Parameters(new Parameter[]
                {
                    new IntParam(MINK, 5) {Help = "Minimum number of peptides per protein required in computing the global/reference standard-deviation in the outlier-filtering step"},
                    new DoubleParam(BFDRCUT, 0.05) {Help = "Statistical significance cut-off." }
                })
        };

        public static Parameter[] GetFiltering()
        {
            return new Parameter[]
            {
                GetOutlierRM()
            };
        }

        public static DoubleParam GetLeftB()
        {
            return new DoubleParam(LEFTB, 0.2) { Help = "Lower bound of the percentile of log ratios from which the null distribution is estimated from." };
        }

        public static DoubleParam GetRightB()
        {
            return new DoubleParam(RIGHTB, 0.8) { Help = "Upper bound of the percentile of log ratios from which the null distribution is estimated up to." };
        }

        public static Parameter[] GetNullDistr()
        {
            return new Parameter[]
            {
                GetLeftB(),
                GetRightB()
            };
        }

        //EBConv from here onwards

        public static string GroupString(int i, string str)
        {
            return "Group " + i.ToString() + " " + str;
        }

        public static string GroupvsGroup(int i, int j)
        {
            return "Group " + i.ToString() + "/Group " + j.ToString();
        }

        public static string[] GetContrastOptions(int option, int n)
        {
            string[] contrast = new string[2 * (n - 1) + 1];
            contrast[2 * (n - 1)] = "None";
            int i = 1;
            int index = 0;
            while (i <= n)
            {
                if (option != i)
                {
                    contrast[index++] = GroupvsGroup(option, i);
                    contrast[index++] = GroupvsGroup(i, option);
                }
                i += 1;
            }
            return contrast;
        }

        public static SingleChoiceWithSubParams GetGroups(IMatrixData mdata)
        {
            string[] choiceList = new string[10];
            //be careful about the difference between Parameter and Parameters!
            Parameters[] subParamList = new Parameters[10];
            int parNameWidth = 100;
            int totWidth = 600;

            for (int i = 0; i < 10; i++)
            {
                choiceList[i] = (i + 1).ToString();
                //2*(i+1)+i+1
                Parameter[] subsubParam = new Parameter[3 * (i + 1)];
                //how many data for each group when equally split
                int split = ArrayUtils.Concat(mdata.ColumnNames, mdata.NumericColumnNames).Length / (i + 1);
                for (int j = 0; j <= i; j++)
                {
                    subsubParam[2 * j] = new StringParam(GroupString(j + 1, "Label"));
                    int[] def = split != 0 ? Enumerable.Range(j * split, split).ToArray() : new int[0];
                    //fill
                    subsubParam[2 * j + 1] = GetData(mdata, GroupString(j + 1, "Data"), "", def);

                    //add contrast here
                    subsubParam[2 * (i + 1) + j] = new SingleChoiceParam(GroupString(j + 1, CONTRAST), 0)
                    {
                        Values = GetContrastOptions(j + 1, i + 1),//new string[] { "1", "2", "3" },//
                        //fill
                        Help = ""
                    };
                }
                subParamList[i] = new Parameters(subsubParam);
                
            }

            return new SingleChoiceWithSubParams(GROUPNUM, 0)
            {
                Values = choiceList,
                SubParams = subParamList,
                ParamNameWidth = parNameWidth,
                TotalWidth = totWidth,
                //fill
                Help = ""
            };
        }

        public static IntParam GetMinSample()
        {
            return new IntParam(MINSAMP, 1) { Help = "Minimum number of samples across each group for the calculation of weighted average intensity. This number should be higher if there are many samples in each group." };
        }

        public static Parameter[] SelectConvData(IMatrixData mdata)
        {
            //fill
            //string name = "Intensity Data";
            //string help = "";

            return new Parameter[]
            {
                GetMinPep(),
                GetMinSample(),
                GetPeptideIDColumn(mdata),
                GetProteinIDColumn(mdata),
                GetGroups(mdata),
            };
        }


        //public static Parameter[] GetConvEBData()
        //{
        //    return new Parameter[]
        //    {

        //        GetExpDesign(),
        //        GetDataForm(),
        //        GetOutlierRM(),
        //        GetLeftB(),
        //        GetRightB()
        //    };
        //}




    }
}
