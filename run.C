// #include "StRoot/StPhiKKMaker/StPhiKKMaker.h"
// #include "StPhiKKMaker/StPhiKKMaker.h"
class StMaker;
class StPhiKKMaker;

void run(){

    gSystem->Load( "St_base" );
    gSystem->Load( "StChain" );
    gSystem->Load( "StEvent" );
    gSystem->Load( "libStDb_Tables.so" );
    gSystem->Load( "libgeometry_Tables.so" );
    gSystem->Load( "StEmcUtil" );


    gSystem->Load( "StStrangeMuDstMaker" );
    gSystem->Load( "StMuDSTMaker" );

    // gSystem->Load( "StMaker" );
    gSystem->Load( "StUpcDst" );
    gSystem->Load( "StPhiKKMaker" );

    StPhiKKMaker* mk = new StPhiKKMaker();
    // StUPCEvent * p = 0;

    mk->Make();

}
