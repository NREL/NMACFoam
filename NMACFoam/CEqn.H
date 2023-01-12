{

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi, C)
      - fvm::laplacian(diffcoeff, C)
    );

    CEqn.relax();

    fvOptions.constrain(CEqn);
    
    fvOptions.correct(C);

    CEqn.solve();

}
