using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text.RegularExpressions;

namespace PluginEBprot
{
    public static class ExternalProcess
    {
        public const int GRPTASK = 0;
        public const int EBPROTTASK = 1;
        public static int RunExe(string exe, string arg, string workingDir, Action<string> status, Action<int> progress, int task, out string errorString)
        {
            errorString = null;
            var externalProcessInfo = new ProcessStartInfo
            {
                FileName = exe,
                Arguments = arg,
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true,
                WorkingDirectory = @workingDir
            };
            var process = new Process { StartInfo = externalProcessInfo };
            var outputData = new List<string>();
            process.OutputDataReceived += (sender, output) =>
            {
                status(output.Data);
                outputData.Add(output.Data);
            };
            var errorData = new List<string>();
            process.ErrorDataReceived += (sender, error) =>
            {
                //Debug.WriteLine(error.Data);
                errorData.Add(error.Data);
            };
            process.Start();
            process.BeginErrorReadLine();
            process.BeginOutputReadLine();
            process.WaitForExit();
            int exitCode = process.ExitCode;
            //Debug.WriteLine($"Process exited with exit code {exitCode}");
            if (exitCode != 0)
            {
                var statusString = String.Join("\n", outputData);
                var errString = String.Join("\n", errorData);
                errorString = String.Concat("Output\n", statusString, "\n", "Error\n", errString);
            }
            process.Dispose();
            if (task == GRPTASK)
            {
                progress(50);
            }
            else if (task == EBPROTTASK)
            {
                progress(100);
            }
            return exitCode;
        }

        public static int RunEBprot(string workingDir, Action<string> status, Action<int> progress, out string errorString)
        {
            string exe = System.IO.Path.Combine(Directory.GetCurrentDirectory(), @".\bin\EBprotInstallations\EBprotV2.exe");
            string args = "\"" + System.IO.Path.Combine(@workingDir, Utils.EBPROTINPUTDATAFILE) + "\"" + " \"" + System.IO.Path.Combine(@workingDir, Utils.EBPROTINPUTPARAMFILE) + "\"";
            //string args = "\"" + System.IO.Path.Combine(@workingDir, @".\weighted_grpcomparisons.txt") + "\"" + " \"" + System.IO.Path.Combine(@workingDir, @".\input_EBprotV2") + "\"";
            return RunExe(exe, args, workingDir, status, progress, EBPROTTASK, out errorString);

        }
        public static int RunGRP(string workingDir, Action<string> status, Action<int> progress, out string errorString)
        {
            string exe = System.IO.Path.Combine(Directory.GetCurrentDirectory(), @".\bin\EBprotInstallations\EBprot.MakeGrpData.exe");
            string args = "\"" + System.IO.Path.Combine(@workingDir, Utils.GRPINPUTDATAFILE) + "\"" + " \"" + System.IO.Path.Combine(@workingDir, Utils.GRPINPUTPARAMFILE) + "\"";
            return RunExe(exe, args, workingDir, status, progress, GRPTASK, out errorString);
        }
    }
}
