class StMaker;
class StPairDstMaker;

void runPairDst() {
    gSystem->Load("St_base");
    gSystem->Load("StChain");
    gSystem->Load("StEvent");
    gSystem->Load("libStDb_Tables.so");
    gSystem->Load("libgeometry_Tables.so");
    gSystem->Load("StEmcUtil");

    gSystem->Load("StStrangeMuDstMaker");
    gSystem->Load("StMuDSTMaker");

    gSystem->Load("StUpcDst");
    gSystem->Load("StPairDstMaker");

    StPairDstMaker* mk = new StPairDstMaker("PairDstMaker");
    mk->SetInputFile("/path/to/your/input/file");
    mk->SetOutputFile("output_femtopair.root");

    mk->Init();
    mk->Make();
    mk->Finish();

    delete mk;
}
