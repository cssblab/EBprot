using PerseusApi.Matrix;
using PerseusApi.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PluginEBprot;
using BaseLibS.Graph.Image;
using BaseLibS.Graph;
using BaseLibS.Param;
using PerseusApi.Document;
using PerseusApi.Generic;

namespace PluginEBprot.EBprot
{
    class MatrixProcessing : PluginEBprot.MatrixProcessing, IMatrixProcessing
    {
        public override string Name => "EBprotV2";

        //fill
        public override string Description => "";

        public override float DisplayRank => 2;

        //fill
        public override Bitmap2 DisplayImage => Utils.GetImage("EBlogo.png");

        public override Parameters GetParameters(IMatrixData mdata, ref string errString)
        {

            Parameters parameters = new Parameters
                (
                    EBParameters.GetWorkingDir()
                );
            parameters.AddParameterGroup(EBParameters.GetAboutEBData(), "About Data", false);
            parameters.AddParameterGroup(EBParameters.GetFiltering(), "Filtering", false);
            parameters.AddParameterGroup(EBParameters.GetNullDistr(), "Null Distribution Estimate", false);
            parameters.AddParameterGroup(EBParameters.SelectEBData(mdata), "Select Data", false);

            return parameters;
        }

        public override void ProcessData(IMatrixData mdata, Parameters param, ref IMatrixData[] supplTables, ref IDocumentData[] documents, ProcessInfo processInfo)
        {
            string workingDirectory = param.GetParam<string>("Working Directory").Value;

            if (Utils.WriteEBInputParam(param, workingDirectory) != 0)
            {
                processInfo.ErrString = "Unable To Process the Given Parameters";
                return;
            }

            if (Utils.WriteInputFiles(mdata, param, workingDirectory, out string errString1, Utils.EBTASK) != 0)
            {
                processInfo.ErrString = errString1;
                return;
            }

            if (ExternalProcess.RunEBprot(workingDirectory, processInfo.Status, processInfo.Progress, out string processInfoErrString) != 0)
            {
                processInfo.ErrString = processInfoErrString;
                return;
            }

            Utils.LoadOutput(mdata, workingDirectory);

            processInfo.Progress(0);
        }
    }
    
}
