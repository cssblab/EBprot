using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BaseLibS.Param;
using PerseusApi.Matrix;
using PerseusApi.Utils;
using BaseLibS.Parse;
using BaseLibS.Num;
using PerseusApi.Generic;
using System.IO;
using System.Text.RegularExpressions;
using BaseLibS.Graph;
using System.Reflection;
using BaseLibS.Graph.Image;

namespace PluginEBprot
{
    public static class Utils
    {
        public const int CONVTASK = 0;
        public const int EBTASK = 1;
        public const string EBPROTINPUTPARAMFILE = "input_EBprotV2";
        public const string GRPINPUTPARAMFILE = "input_MakeGroupData";
        public const string EBPROTINPUTDATAFILE = "weighted_grpcomparisons.txt";
        public const string GRPINPUTDATAFILE = "data_input.txt";
        public const string OUTPUTFILE = "EBprot_results.txt";

        //index of bases
        public const int RAW = 0;
        public const int LOG2 = 1;
        public const int LN = 2;
        public const int LOG10 = 3;
        public const int CUSTOM = 4;

        //-1 to signal it's untransformed
        public const double RAWBASE = -1;
        //0 means no need to change
        public const double LOG2BASE = 0;
        public const double LNBASE = Math.E;
        public const double LOG10BASE = 10;
        public const double ERRBASE = -2;

        public static Bitmap2 GetImage(string file)
        {
            //remember to change the image files Build Action to Embedded Resource
            Assembly thisExe = Assembly.GetExecutingAssembly();
            //path is default namespace + folder, where each is separated by '.'
            Stream file1 = thisExe.GetManifestResourceStream("PluginEBprot.Resources." + file);
            if (file1 == null)
            {
                return null;
            }
            Bitmap2 bm = Image2.ReadImage(file1);
            file1.Close();
            return bm;
        }

        public static double GetBase(ParameterWithSubParams<int> form)
        {
            
            int formType = form.Value;

            if (formType == RAW) return RAWBASE;

            if (formType == LN) return LOG2BASE;

            if (formType == LOG10)
            {
                return LOG10BASE;
            }
            else if (formType == LOG2)
            {
                return LOG2BASE;
            }
            else if (formType == CUSTOM)
            {
                Parameters subParams = form.GetSubParameters();
                return subParams.GetParam<double>("Base").Value;
            }
            else
            {
                return ERRBASE;
            }

        }

        public static string GetDesign(int val)
        {
            switch(val)
            {
                case 0:
                    return "independent";
                case 1:
                    return "replicate";
                case 2:
                    return "timecourse";
                default:
                    return null;
            }
        }

        public static string[] GetEBInputParamLines(Parameters param, string groupLabels = null, List<int> constrastIndices = null)
        {
            //shared by both EBprot and Conv Tasks
            string design = GetDesign(param.GetParam<int>(EBParameters.EXPDESIGN).Value);
            string needLog2Transformed = "";
            string labelString = "";
            //EBprot when contrastString is null
            if (constrastIndices ==null)
            {
                needLog2Transformed = GetBase(param.GetParamWithSubParams<int>(EBParameters.DATAFORM)) == 0 ? "false" : "true";
                ParameterWithSubParams<int> labels = param.GetParamWithSubParams<int>(EBParameters.LABELNUM);
                //labels.Value give the index from 1 to 10, so add 1
                int labNum = labels.Value + 1;
                for (int i = 0; i < labNum; i++)
                {
                    labelString += (labels.GetSubParameters().GetParam<string>(EBParameters.LabelString(i + 1, "Name")).Value + " ");
                }
            }
            else
            {
                //already transformed
                needLog2Transformed = "false";
                string[] groupLab = groupLabels.Split(' ');
                //there's an extra space
                int groupNum = groupLab.Length - 1;
                for (int i = 0; i < constrastIndices.Count; i++)
                {
                    int groupIdx1 = (constrastIndices[i]) / groupNum;
                    int groupIdx2 = (constrastIndices[i]) % groupNum;
                    labelString += (groupLab[groupIdx1] + "vs" + groupLab[groupIdx2] + " ");
                }
            }
            
            //shared
            ParameterWithSubParams<bool> outlierFilt = param.GetParamWithSubParams<bool>(EBParameters.OUTFILT);
            string outlierRM = outlierFilt.Value ? "true" : "false";
            int minPep = param.GetParam<int>(EBParameters.MINPEP).Value;
            double leftB = param.GetParam<double>(EBParameters.LEFTB).Value;
            double rightB = param.GetParam<double>(EBParameters.RIGHTB).Value;
            double BFDR = param.GetParam<double>(EBParameters.BFDRCUT).Value;

            string[] lines =
            {
                string.Format("EXPERIMENTAL_DESIGN = {0}", design),
                string.Format("LOG2_TRANSFORM = {0}", needLog2Transformed),
                string.Format("OUTLIER_RM = {0}", outlierRM),
                //placeholders for outfilt true
                //indices 3
                "",
                string.Format("BFDR_THRESHOLD = {0}", BFDR),
                string.Format("LABELS = {0}",labelString),
                string.Format("MIN_PEP = {0}", minPep),
                string.Format("LEFT_B = {0}", leftB),
                string.Format("RIGHT_B = {0}", rightB)
            };
            if (outlierFilt.Value)
            {
                int minK = outlierFilt.GetSubParameters().GetParam<int>(EBParameters.MINK).Value;
                //double BFDR = outlierFilt.GetSubParameters().GetParam<double>(EBParameters.BFDRCUT).Value;
                lines[3] = string.Format("MIN_K = {0}", minK);
                //lines[4] = string.Format("BFDR_THRESHOLD = {0}", BFDR);
            }
            return lines;
        }

        public static int WriteEBInputParam(Parameters param, string workingDir)
        {
            string[] lines = GetEBInputParamLines(param);
            string inputLocation = Path.Combine(workingDir, EBPROTINPUTPARAMFILE);
            File.WriteAllLines(inputLocation, lines);
            return 0;
        }

        //modified from EBParameters.GetContrastOptions
        //returns index for the contrast matrix
        //option refers to the group number
        public static int GetContrastIdx(int option, int n, int tarIdx)
        {
            int i = 1;
            int index = 0;
            
            //to indicate none
            if (tarIdx == 2 * (n - 1)) return -1;
            while (i <= n)
            {
                if (option != i)
                {
                    //row option columnn i
                    if (index++ == tarIdx) return (option - 1) * n + (i - 1);
                    //row i column option
                    if (index++ == tarIdx) return (i - 1) * n + (option - 1);

                }
                i += 1;
            }
            //did not find a matching targetIdx so error
            return -2;
        }
        
        public static string GetContrastMatrixString(Parameters param, int groupNum, out List<int> indices)
        {
            indices = new List<int>();
            int gridNum = groupNum * groupNum;
            string[] contrastMtx = new string[gridNum * 2];
            for (int i = 0; i < groupNum; i++)
            {
                int idx = param.GetParam<int>(EBParameters.GroupString(i + 1, EBParameters.CONTRAST)).Value;
                int contrastIdx = GetContrastIdx(i + 1, groupNum, idx);
                if (contrastIdx == -2)
                {
                    //error
                    indices = null;
                    return null;
                }
                else if (contrastIdx != -1)
                {
                    indices.Add(contrastIdx);
                    contrastMtx[2 * contrastIdx] = "1";
                }
            }
            for (int i = 0; i < groupNum; i++)
            {
                for (int j = 0; j < groupNum; j++)
                {
                    int idx = (groupNum * i + j) * 2;
                    if (i == j) contrastMtx[idx] = "-";
                    else
                    {
                        if (contrastMtx[idx] == null) contrastMtx[idx] = "0";
                    }
                }
            }
            //filling out space and newline characters
            for (int i = 0; i < gridNum; i++)
            {
                contrastMtx[2 * i + 1] = ((i + 1) % groupNum == 0) ? "\n" : " ";
            }
            indices.Sort();
            return String.Join("",contrastMtx);
        }

        public static string[] GetGrpInputParamLines(Parameters param, out string groupLabels, out List<int> indices)
        {
            ParameterWithSubParams<int> groups = param.GetParamWithSubParams<int>(EBParameters.GROUPNUM);
            //groups.Value give the index from 1 to 10, so add 1
            int groupNum = groups.Value + 1;
            string groupSizes = "";
            groupLabels = "";
            for (int i = 0; i < groupNum; i++)
            {
                int size = groups.GetSubParameters().GetParam<int[]>(EBParameters.GroupString(i + 1, "Data")).Value.Length;
                groupSizes += (size.ToString() + " ");
                string label = groups.GetSubParameters().GetParam<string>(EBParameters.GroupString(i + 1, "Name")).Value;
                groupLabels += (label + " ");
            }
            int minSamp = param.GetParam<int>(EBParameters.MINSAMP).Value;
            string isLog2 = GetBase(param.GetParamWithSubParams<int>(EBParameters.DATAFORM)) == 0 ? "true" : "false";
            string contrastString = GetContrastMatrixString(groups.GetSubParameters(), groupNum, out indices);

            string[] lines =
            {
                string.Format("GROUP_SIZE = {0}", groupSizes),
                string.Format("GROUP_LABELS = {0}", groupLabels),
                string.Format("MIN_SAMPLE = {0}", minSamp),
                string.Format("IS_LOG = {0}",isLog2),
                string.Format("CONTRAST_MATRIX = \n{0}", contrastString),
            };

            return lines;
        }


        public static int WriteGrpEBInputParam(Parameters param, string workingDir)
        {
            string[] GRPlines = GetGrpInputParamLines(param, out string groupLabels, out List<int> contrastIndices);
            if (contrastIndices == null) return -1;
            string[] EBlines = GetEBInputParamLines(param, groupLabels, contrastIndices);
            string inputGrpLocation = Path.Combine(workingDir, GRPINPUTPARAMFILE);
            string inputEBLocation = Path.Combine(workingDir, EBPROTINPUTPARAMFILE);
            File.WriteAllLines(inputGrpLocation, GRPlines);
            File.WriteAllLines(inputEBLocation, EBlines);
            return 0;
        }

        public static int WriteInputFiles(IMatrixData mdata, Parameters param, string workingDirectory, out string errString, int task = EBTASK)
        {
            string inputFile = task == EBTASK ? EBPROTINPUTDATAFILE : GRPINPUTDATAFILE;
            string tempFile = Path.Combine(workingDirectory, ("Temp" + inputFile));
            string finFile = Path.Combine(workingDirectory, inputFile);

            IMatrixData mCopy = (IMatrixData)mdata.Clone();

            int[] nameInd = new[] { param.GetParam<int>(EBParameters.PEPID).Value, param.GetParam<int>(EBParameters.PROTID).Value };
            List<int> valIdx = new List<int>();

            if (task == EBTASK)
            {
                ParameterWithSubParams<int> labels = param.GetParamWithSubParams<int>(EBParameters.LABELNUM);
                int num = labels.Value + 1;
                for (int i = 0; i < num; i++)
                {
                    valIdx.AddRange(labels.GetSubParameters().GetParam<int[]>(EBParameters.LabelString(i + 1, "Data")).Value);
                }
            }
            else if (task == CONVTASK)
            {
                ParameterWithSubParams<int> groups = param.GetParamWithSubParams<int>(EBParameters.GROUPNUM);
                int num = groups.Value + 1;
                for (int i = 0; i < num; i++)
                {
                    valIdx.AddRange(groups.GetSubParameters().GetParam<int[]>(EBParameters.GroupString(i + 1, "Data")).Value);
                }
            }
            else
            {
                errString = "Wrong Task #";
                return -1;
            }

            int[] valInd = valIdx.ToArray();
            double baseVal = GetBase(param.GetParamWithSubParams<int>(EBParameters.DATAFORM));

            SetupMDataForInput(mCopy, valInd, nameInd, baseVal);

            try
            {
                PerseusUtils.WriteMatrixToFile(mCopy, tempFile, false);
                RemoveCommentsFromFile(tempFile, finFile);
            }
            catch (Exception e)
            {
                errString = e.ToString();
                return -1;
            }
            errString = null;
            return 0;
        }

        public static void LoadOutput(IMatrixData mdata, Parameters param, string workingDir)
        {
            string filename = Path.Combine(workingDir, OUTPUTFILE);
            char separator = '\t';
            string[] colNames = TabSep.GetColumnNames(filename, 0, PerseusUtils.commentPrefix,
                PerseusUtils.commentPrefixExceptions, null, separator);
            List<string> newColNames = new List<string>(colNames);

            if (!param.GetParam<bool>(EBParameters.MEDRATIO).Value)
            {
                newColNames.RemoveAll(str => str.Contains("MedianlogRatio"));
            }
            if (!param.GetParam<bool>(EBParameters.NUMPEP).Value)
            {
                newColNames.RemoveAll(str => str.Contains("NumPep_"));
            }
            if (!param.GetParam<bool>(EBParameters.PEPREM).Value)
            {
                newColNames.RemoveAll(str => str.Contains("NumPepRm"));
            }
            if (!param.GetParam<bool>(EBParameters.PPSCORE).Value)
            {
                newColNames.RemoveAll(str => str.Contains("PPscore"));
            }
            if (!param.GetParam<bool>(EBParameters.POSTODD).Value)
            {
                newColNames.RemoveAll(str => str.Contains("PostOdds"));
            }

            string[][] cols = TabSep.GetColumns(newColNames.ToArray(), filename, 0, PerseusUtils.commentPrefix,
                PerseusUtils.commentPrefixExceptions, separator);
            int nrows = TabSep.GetRowCount(filename);


            mdata.Clear();
            mdata.Name = "EBprot Result";
            mdata.Values.Init(nrows, 0);
            mdata.SetAnnotationColumns(new List<string>(newColNames), new List<string>(colNames), new List<string[]>(cols), new List<string>(),
                new List<string>(), new List<string[][]>(), new List<string>(), new List<string>(), new List<double[]>(),
                new List<string>(), new List<string>(), new List<double[][]>());

            //no numeric expression list at this point
            int[] numericList = Enumerable.Range(1, newColNames.Count - 1).ToArray();
            StringToNumerical(numericList, mdata);

        }



        //other data handling

        //function obtained from PerseusPluginLib/Rearrange/ChangeColumnType.cs
        private static void StringToNumerical(IList<int> colInds, IMatrixData mdata)
        {
            int[] inds = ArrayUtils.Complement(colInds, mdata.StringColumnCount);
            string[] name = ArrayUtils.SubArray(mdata.StringColumnNames, colInds);
            string[] description = ArrayUtils.SubArray(mdata.StringColumnDescriptions, colInds);
            string[][] str = ArrayUtils.SubArray(mdata.StringColumns, colInds);
            var newNum = new double[str.Length][];
            for (int j = 0; j < str.Length; j++)
            {
                newNum[j] = new double[str[j].Length];
                for (int i = 0; i < newNum[j].Length; i++)
                {
                    if (str[j][i] == null || str[j][i].Length == 0)
                    {
                        newNum[j][i] = double.NaN;
                    }
                    else
                    {
                        string x = str[j][i];
                        double d;
                        bool success = double.TryParse(x, out d);
                        newNum[j][i] = success ? d : double.NaN;
                    }
                }
            }
            mdata.NumericColumnNames.AddRange(name);
            mdata.NumericColumnDescriptions.AddRange(description);
            mdata.NumericColumns.AddRange(newNum);
            mdata.StringColumns = ArrayUtils.SubList(mdata.StringColumns, inds);
            mdata.StringColumnNames = ArrayUtils.SubList(mdata.StringColumnNames, inds);
            mdata.StringColumnDescriptions = ArrayUtils.SubList(mdata.StringColumnDescriptions, inds);
        }

        //function obtained from PerseusPluginLib/Rearrange/ChangeColumnType.cs
        private static void NumericToString(IList<int> colInds, IMatrixData mdata)
        {
            int[] inds = ArrayUtils.Complement(colInds, mdata.NumericColumnCount);
            string[] name = ArrayUtils.SubArray(mdata.NumericColumnNames, colInds);
            string[] description = ArrayUtils.SubArray(mdata.NumericColumnDescriptions, colInds);
            double[][] num = ArrayUtils.SubArray(mdata.NumericColumns, colInds);
            var newString = new string[num.Length][];
            for (int j = 0; j < num.Length; j++)
            {
                newString[j] = new string[num[j].Length];
                for (int i = 0; i < newString[j].Length; i++)
                {
                    double x = num[j][i];
                    newString[j][i] = "" + x;
                }
            }
            mdata.StringColumnNames.AddRange(name);
            mdata.StringColumnDescriptions.AddRange(description);
            mdata.StringColumns.AddRange(newString);
            mdata.NumericColumns = ArrayUtils.SubList(mdata.NumericColumns, inds);
            mdata.NumericColumnNames = ArrayUtils.SubList(mdata.NumericColumnNames, inds);
            mdata.NumericColumnDescriptions = ArrayUtils.SubList(mdata.NumericColumnDescriptions, inds);
        }
        //function obtained from PerseusPluginLib/Rearrange/ChangeColumnType.cs
        private static void ExpressionToNumeric(IList<int> colInds, IMatrixData mdata)
        {
            int[] remainingInds = ArrayUtils.Complement(colInds, mdata.ColumnCount);
            foreach (int colInd in colInds)
            {
                double[] d = ArrayUtils.ToDoubles(mdata.Values.GetColumn(colInd));
                mdata.AddNumericColumn(mdata.ColumnNames[colInd], mdata.ColumnDescriptions[colInd], d);
            }
            mdata.ExtractColumns(remainingInds);
        }
        //function obtained from PerseusPluginLib/Rearrange/ChangeColumnType.cs
        private static void StringToExpression(IList<int> colInds, IMatrixData mdata)
        {
            int[] inds = ArrayUtils.Complement(colInds, mdata.StringColumnCount);
            string[] names = ArrayUtils.SubArray(mdata.StringColumnNames, colInds);
            string[] descriptions = ArrayUtils.SubArray(mdata.StringColumnDescriptions, colInds);
            string[][] str = ArrayUtils.SubArray(mdata.StringColumns, colInds);
            var newEx = new float[str.Length][];
            for (int j = 0; j < str.Length; j++)
            {
                newEx[j] = new float[str[j].Length];
                for (int i = 0; i < newEx[j].Length; i++)
                {
                    float f;
                    bool success = float.TryParse(str[j][i], out f);
                    newEx[j][i] = success ? f : float.NaN;
                }
            }
            float[,] newExp = new float[mdata.RowCount, mdata.ColumnCount + str.Length];
            float[,] newQual = new float[mdata.RowCount, mdata.ColumnCount + str.Length];
            bool[,] newIsImputed = new bool[mdata.RowCount, mdata.ColumnCount + str.Length];
            for (int i = 0; i < mdata.RowCount; i++)
            {
                for (int j = 0; j < mdata.ColumnCount; j++)
                {
                    newExp[i, j] = mdata.Values.Get(i, j);
                    newQual[i, j] = mdata.Quality.Get(i, j);
                    newIsImputed[i, j] = mdata.IsImputed[i, j];
                }
                for (int j = 0; j < newEx.Length; j++)
                {
                    newExp[i, j + mdata.ColumnCount] = newEx[j][i];
                    newQual[i, j + mdata.ColumnCount] = float.NaN;
                    newIsImputed[i, j + mdata.ColumnCount] = false;
                }
            }
            mdata.Values.Set(newExp);
            mdata.Quality.Set(newQual);
            mdata.IsImputed.Set(newIsImputed);
            mdata.ColumnNames.AddRange(names);
            mdata.ColumnDescriptions.AddRange(descriptions);
            mdata.StringColumns = ArrayUtils.SubList(mdata.StringColumns, inds);
            mdata.StringColumnNames = ArrayUtils.SubList(mdata.StringColumnNames, inds);
            mdata.StringColumnDescriptions = ArrayUtils.SubList(mdata.StringColumnDescriptions, inds);
            for (int i = 0; i < mdata.CategoryRowCount; i++)
            {
                mdata.SetCategoryRowAt(ExtendCategoryRow(mdata.GetCategoryRowAt(i), str.Length), i);
            }
            for (int i = 0; i < mdata.NumericRows.Count; i++)
            {
                mdata.NumericRows[i] = ExtendNumericRow(mdata.NumericRows[i], str.Length);
            }
        }
        //function obtained from PerseusPluginLib/Rearrange/ChangeColumnType.cs
        private static double[] ExtendNumericRow(IList<double> numericRow, int add)
        {
            var result = new double[numericRow.Count + add];
            for (int i = 0; i < numericRow.Count; i++)
            {
                result[i] = numericRow[i];
            }
            for (int i = numericRow.Count; i < numericRow.Count + add; i++)
            {
                result[i] = double.NaN;
            }
            return result;
        }
        //function obtained from PerseusPluginLib/Rearrange/ChangeColumnType.cs
        private static string[][] ExtendCategoryRow(IList<string[]> categoryRow, int add)
        {
            var result = new string[categoryRow.Count + add][];
            for (int i = 0; i < categoryRow.Count; i++)
            {
                result[i] = categoryRow[i];
            }
            for (int i = categoryRow.Count; i < categoryRow.Count + add; i++)
            {
                result[i] = new string[0];
            }
            return result;
        }

        public static void RemoveCommentsFromFile(string filename, string destFile)//, int lines_to_delete)
        {
            using (StreamReader reader = new StreamReader(filename))
            using (StreamWriter writer = new StreamWriter(destFile))
            {
                writer.WriteLine(reader.ReadLine());
                string line = reader.ReadLine();
                while (line.Contains("#!"))
                {
                    line = reader.ReadLine();
                }
                writer.WriteLine(line);
                while ((line = reader.ReadLine()) != null)
                    writer.WriteLine(line);
            }
            File.Delete(filename);
        }


        public static void SetupMDataForInput(IMatrixData data, int[] columnIndx, int[] nameInd, double baseVal)
        {

            data.StringColumns = ArrayUtils.SubList(data.StringColumns, nameInd);
            data.StringColumnNames = ArrayUtils.SubList(data.StringColumnNames, nameInd);
            data.StringColumnDescriptions = ArrayUtils.SubList(data.StringColumnDescriptions, nameInd);

            List<int> toConvert = new List<int>();
            List<int> numList = new List<int>();
            int expressInd = 0;
            foreach (int i in columnIndx)
            {
                if (i < data.ColumnCount)
                {
                    toConvert.Add(i);
                    numList.Add(data.NumericColumnCount + expressInd);
                    expressInd += 1;
                }
                else
                {
                    numList.Add(i - data.ColumnCount);
                }
            }

            int[] numArr = numList.ToArray();
            //convert expression to numeric
            data.ExtractColumns(toConvert.ToArray());
            ExpressionToNumeric(Enumerable.Range(0, data.ColumnCount).ToArray(), data);

            data.NumericColumns = ArrayUtils.SubList(data.NumericColumns, numArr);

            //change data form depending whether needed
            if (baseVal > 0)
            {
                foreach (int col in numArr)
                {
                    for (int i = 0; i < data.RowCount; i++)
                    {
                        data.NumericColumns[col][i] = Math.Pow(baseVal, data.NumericColumns[col][i]);
                    }
                }
            }

            data.NumericColumnNames = ArrayUtils.SubList(data.NumericColumnNames, numArr);
            data.NumericColumnDescriptions = ArrayUtils.SubList(data.NumericColumnDescriptions, numArr);

            NumericToString(Enumerable.Range(0, numArr.Length).ToArray(), data);


            for (int j = 0; j < data.StringColumnCount; j++)
            {
                for (int i = 0; i < data.RowCount; i++)
                {
                    data.StringColumns[j][i] = string.Equals(data.StringColumns[j][i], "NaN") ? "NA" : data.StringColumns[j][i];
                }
            }

            //clearing irrelevant info
            data.ClearMultiNumericColumns();
            data.ClearMultiNumericRows();
            data.ClearCategoryColumns();
            data.ClearCategoryRows();

        }

    }
}
