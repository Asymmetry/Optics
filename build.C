void build()
{
    gROOT->ProcessLine(".L LOpticsOpt.C+");
    gROOT->ProcessLine(".L ROpticsOpt.C+");
}
