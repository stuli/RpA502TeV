
void WriteToLog() {


  ofstream logFile;
  logFile.open("log.txt", fstream::in | fstream::out | fstream::app);

  logFile << "**UNABLE TO GET A SATISFACTORY FIT IN THIS BIN!**" << endl;

}
