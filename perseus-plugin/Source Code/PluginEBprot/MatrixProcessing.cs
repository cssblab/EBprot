using BaseLibS.Graph;
using BaseLibS.Num;
using BaseLibS.Param;
using BaseLibS.Parse;
using PerseusApi.Document;
using PerseusApi.Generic;
using PerseusApi.Matrix;
using PerseusApi.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace PluginEBprot
{
    public abstract class MatrixProcessing : IMatrixProcessing
    {
        public abstract string Name { get; }

        public abstract string Description { get; }

        public virtual float DisplayRank => 1;

        public virtual bool IsActive => true;

        public virtual int GetMaxThreads(global::BaseLibS.Param.Parameters parameters) => 1;

        public virtual bool HasButton => false;

        public virtual Bitmap2 DisplayImage => null;

        //insert github plugin in the future
        public virtual string Url => "https://www.ncbi.nlm.nih.gov/pubmed/25913743";

        public virtual string Heading => "EBprot";

        public virtual string HelpOutput { get; }

        public virtual string[] HelpSupplTables { get; }

        public virtual int NumSupplTables { get; }

        public virtual string[] HelpDocuments { get; }

        public virtual int NumDocuments { get; }

        public virtual Parameters GetParameters(IMatrixData mdata, ref string errString)
        {
            Parameters parameters = new Parameters();
            return parameters;
        }
        public virtual void ProcessData(IMatrixData mdata, global::BaseLibS.Param.Parameters param, ref IMatrixData[] supplTables, ref IDocumentData[] documents, ProcessInfo processInfo)
        {

        }
    }
}
