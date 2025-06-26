// void runPairDst() {
//     gSystem->Load("St_base");
//     gSystem->Load("StChain");
//     gSystem->Load("StEvent");
//     gSystem->Load("libStDb_Tables.so");
//     gSystem->Load("libgeometry_Tables.so");
//     gSystem->Load("StEmcUtil");

//     gSystem->Load("StStrangeMuDstMaker");
//     gSystem->Load("StMuDSTMaker");

//     gSystem->Load("StUpcDst");
//     gSystem->Load("StPairDstMaker");

//     StPairDstMaker* mk = new StPairDstMaker("PairDstMaker");
//     mk->SetInputFileList("mid14.lis");
//     mk->SetOutputFile("output_femtopair.root");
//     std::vector<Int_t> triggers;
//     triggers.push_back(450701);
//     triggers.push_back(450711);
//     mk->SetTriggerIds( triggers );  // Set the desired trigger ID

//     mk->Init();
//     mk->Make();
//     mk->Finish();

//     delete mk;
// }
//
#include <TString.h>
#include <TSystem.h>
#include <vector>

void runPairDst(const char* inputFileList = "mid14.lis", const char* outputFileName = "output_femtopair.root") {
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
    mk->SetInputFileList(inputFileList);
    mk->SetOutputFile(outputFileName);
    std::vector<Int_t> triggers;

    //run14 triggers
    // triggers.push_back(450701); //upc main
    // triggers.push_back(450711);

    // triggers.push_back(450707); //upc main p
    // triggers.push_back(450717);
    // triggers.push_back(450727);

    // // run10 triggers
    // triggers.push_back(260750);

    // run 11 triggers
    // triggers.push_back(4); //upc main
    // triggers.push_back(350007); //upc main
    // triggers.push_back(350017); //upc main
    //
    //run19 triggers
    // triggers.push_back(700001); //zdc min bias

    //run21 triggers
    triggers.push_back(860702); //upc inclusive
    // triggers.push_back(860701); //upc single zdc

    mk->SetTriggerIds(triggers);  // Set the desired trigger ID

    mk->Init();
    mk->Make();
    mk->Finish();

    delete mk;
}

// This function will be called to run the macro with command line arguments
void runPairDstWrapper(const char* inputFileList, const char* outputFileName) {
    runPairDst(inputFileList, outputFileName);
}
