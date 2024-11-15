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
    mk->SetInputFileList("mid14.lis");
    mk->SetOutputFile("output_femtopair.root");
    std::vector<Int_t> triggers = {450701, 450711};
    mk->SetTriggerIds( triggers );  // Set the desired trigger ID

    mk->Init();
    mk->Make();
    mk->Finish();

    delete mk;
}
