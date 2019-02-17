/***************************************************************************
 *
 * Authors:     Jeison Méndez García (jmendez@utp.edu.co)
 *
 * Instituto de Investigaciones en Matemáticas Aplicadas y en Sistemas -- IIMAS
 * Universidad Nacional Autónoma de México -UNAM
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "angular_assignment_mag.h"

void ProgAngularAssignmentMag::defineParams()
{
    //usage
    addUsageLine("Generates a list of candidates for angular assignment for each experimental image");
    //params
    addParamsLine("   -i <md_file>               : Metadata file with input experimental projections");
    addParamsLine("   -ref <md_file>             : Metadata file with input reference projections");
    addParamsLine("  [-odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [-sym <symfile=c1>]         : Enforce symmetry in projections");
}

// Read arguments ==========================================================
void ProgAngularAssignmentMag::readParams()
{
    fnIn = getParam("-i");
    fnRef = getParam("-ref");
    fnDir = getParam("-odir");
    fnSym = getParam("-sym");
}

// Show ====================================================================
void ProgAngularAssignmentMag::show()
{
    if (verbose > 0)
    {
        std::cout << "Input metadata              : "  << fnIn        << std::endl;
        std::cout << "Input references            : "  << fnRef       << std::endl;
        std::cout << "Output directory            : "  << fnDir       << std::endl;
        if (fnSym != "")
            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
    }
}

// Main routine ------------------------------------------------------------
void ProgAngularAssignmentMag::run()
{
    produceSideinfo(); // read metaData file

    FileName fnImgExp;
    FileName fnImgRef;
    MDRow rowExp, rowRef;
    int countInImg = 0, countRefImg = 0;

    // reading input stack image
    int sizeMdIn = mdIn.size();
    size_t Zdim, Ndim;
    getImageSize(mdIn,Xdim,Ydim,Zdim,Ndim);

    // reading reference stack image

    int sizeMdRef = mdRef.size();

    // passing images to Image and then to MultidimArray structure
    const size_t n_bands = 16;
    const size_t startBand = 5;
    const size_t finalBand = n_bands + startBand;
    const size_t n_rad = size_t(Xdim/2 + 0.5);
    size_t n_ang = size_t(180);
    size_t n_ang2 = 2*n_ang;


    // experimental image related
    Image<double>                           ImgIn;
    MultidimArray<double>                   MDaIn(Ydim,Xdim);
    MultidimArray< std::complex<double> >   MDaInF  ;
    MultidimArray< std::complex<double> >   MDaInF2 ;
    MultidimArray<double>                   MDaInFM ;
    MultidimArray<double>                   MDaInFMs;
    MultidimArray<double>                   MDaInFMs_polarPart(n_bands, n_ang2);
    MultidimArray< std::complex<double> >   MDaInFMs_polarF;

    // reference image related
    Image<double>                           ImgRef;
    MultidimArray<double>                   MDaRef(Ydim,Xdim);
    MultidimArray< std::complex<double> >   MDaRefF  ;
    MultidimArray< std::complex<double> >   MDaRefF2 ;
    MultidimArray<double>                   MDaRefFM ;
    MultidimArray<double>                   MDaRefFMs;
    MultidimArray<double>                   MDaRefFMs_polarPart(n_bands, n_ang2);
    MultidimArray< std::complex<double> >   MDaRefFMs_polarF;

    // try to storage all data related to reference images in memory
    for (int k = 0; k < sizeMdRef; k++){
        // reading image
        mdRef.getRow(rowRef, size_t(k+1) );
        rowRef.getValue(MDL_IMAGE, fnImgRef);
        // processing reference image
        ImgRef.read(fnImgRef);
        MDaRef = ImgRef();
        vecMDaRef.push_back(MDaRef);
        _applyFourierImage2(MDaRef, MDaRefF);
        vecMDaRefF.push_back(MDaRefF);
        transformerImage.getCompleteFourier(MDaRefF2);
        _getComplexMagnitude(MDaRefF2, MDaRefFM);
        completeFourierShift(MDaRefFM, MDaRefFMs);
        MDaRefFMs_polarPart = imToPolar(MDaRefFMs,startBand,finalBand,n_bands, n_rad, n_ang2);
        _applyFourierImage2(MDaRefFMs_polarPart, MDaRefFMs_polarF, n_ang);
        vecMDaRefFMs_polarF.push_back(MDaRefFMs_polarF);

    }

    _delayAxes(Ydim, Xdim, n_ang);

    // CCV result matrix
    MultidimArray<double>                   ccMatrixRot;
    MultidimArray<double>                   ccVectorRot;
    std::vector<double>                     cand; // rotation candidates
    int                                     peaksFound = 0; // peaksFound in ccVectorRot
    double                                  tempCoeff;

    // candidates for each loop
    std::vector<unsigned int>               candidatesFirstLoop(sizeMdRef,0);
    std::vector<unsigned int>               Idx(sizeMdRef,0);
    std::vector<double>                     candidatesFirstLoopCoeff(sizeMdRef,0);
    std::vector<double>                     bestTx(sizeMdRef,0);
    std::vector<double>                     bestTy(sizeMdRef,0);
    std::vector<double>                     bestRot(sizeMdRef,0);

    // main loop, input stack

    size_t idxOut; // index for metadata output file
    mdOut.setComment("experiment for metadata output containing data for reconstruction");

    for (countInImg = 0; countInImg < sizeMdIn; countInImg++ ){
        // read experimental image
        mdIn.getRow(rowExp, size_t(countInImg+1) );
        rowExp.getValue(MDL_IMAGE, fnImgExp);
        printf("\r%d of %d", countInImg+1, sizeMdIn);
        fflush(stdout);
        // processing input image
        ImgIn.read(fnImgExp);
        MDaIn = ImgIn();
        _applyFourierImage2(MDaIn, MDaInF);
        transformerImage.getCompleteFourier(MDaInF2);
        _getComplexMagnitude(MDaInF2, MDaInFM);
        completeFourierShift(MDaInFM, MDaInFMs);
        MDaInFMs_polarPart = imToPolar(MDaInFMs,startBand,finalBand,n_bands, n_rad, n_ang2);
        _applyFourierImage2(MDaInFMs_polarPart, MDaInFMs_polarF, n_ang);

        tempCoeff = -10.0;
        int k = 0;
        double bestCandVar, bestCoeff, Tx, Ty;
        // loop over reference stack
        for(countRefImg = 0; countRefImg < sizeMdRef; countRefImg++){
            // computing relative rotation and traslation
            ccMatrix(MDaInFMs_polarF, vecMDaRefFMs_polarF[countRefImg], ccMatrixRot);
            maxByColumn(ccMatrixRot, ccVectorRot, YSIZE(ccMatrixRot), XSIZE(ccMatrixRot));
            peaksFound = 0;
            std::vector<double>().swap(cand);
            rotCandidates3(ccVectorRot, cand, XSIZE(ccMatrixRot), &peaksFound); 
            bestCand(MDaIn, MDaInF, vecMDaRef[countRefImg], cand, peaksFound, &bestCandVar, &Tx, &Ty, &bestCoeff);
            // all the results are storaged for posterior partial_sort
            Idx[countRefImg] = k++;
            candidatesFirstLoop[countRefImg] = countRefImg+1;
            candidatesFirstLoopCoeff[countRefImg] = bestCoeff;
            bestTx[countRefImg] = Tx;
            bestTy[countRefImg] = Ty;
            bestRot[countRefImg] = bestCandVar;

        }
        // choose nCand of the candidates with best corrCoeff
        int nCand = 25;
        std::partial_sort(Idx.begin(), Idx.begin()+nCand, Idx.end(),
                          [&](int i, int j){return candidatesFirstLoopCoeff[i] > candidatesFirstLoopCoeff[j]; });

        // second loop applies search only over better candidates
        k = 0;
        // candidates second loop
        std::vector<unsigned int>               candidatesFirstLoop2(nCand,0);
        std::vector<unsigned int>               Idx2(nCand,0);
        std::vector<double>                     candidatesFirstLoopCoeff2(nCand,0);
        std::vector<double>                     bestTx2(nCand,0);
        std::vector<double>                     bestTy2(nCand,0);
        std::vector<double>                     bestRot2(nCand,0);
        MultidimArray<double>                   MDaRefTrans;
        for (int i = 0; i < nCand; i++){
            // apply transform to reference images and recompute rotational and traslational parameters
            double rotVal = bestRot[ Idx[i] ];
            double trasXval = bestTx[ Idx[i] ];
            double trasYval = bestTy[ Idx[i] ];
            _applyRotationAndShift(vecMDaRef[ Idx[i] ], rotVal, trasXval, trasYval, MDaRefTrans);
            _applyFourierImage2(MDaRefTrans, MDaRefF);
            transformerImage.getCompleteFourier(MDaRefF2);
            _getComplexMagnitude(MDaRefF2, MDaRefFM);
            completeFourierShift(MDaRefFM, MDaRefFMs);
            MDaRefFMs_polarPart = imToPolar(MDaRefFMs,startBand,finalBand,n_bands, n_rad, n_ang2);
            _applyFourierImage2(MDaRefFMs_polarPart, MDaRefFMs_polarF, n_ang);
            // computing relative rotation and traslation
            ccMatrix(MDaInFMs_polarF, MDaRefFMs_polarF, ccMatrixRot);
            maxByColumn(ccMatrixRot, ccVectorRot, YSIZE(ccMatrixRot), XSIZE(ccMatrixRot));
            peaksFound = 0;
            std::vector<double>().swap(cand);
            rotCandidates3(ccVectorRot, cand, XSIZE(ccMatrixRot), &peaksFound);
            bestCand2(MDaIn, MDaInF, MDaRefTrans, cand, peaksFound, &bestCandVar, &Tx, &Ty, &bestCoeff);
            Idx2[i] = k++;
            candidatesFirstLoop2[i] = candidatesFirstLoop[ Idx[i] ];
            candidatesFirstLoopCoeff2[i] = bestCoeff;
            bestTx2[i] = Tx + trasXval;
            bestTy2[i] = Ty + trasYval;
            bestRot2[i] = bestCandVar + rotVal;
        }

        // choose nCand of the candidates with best corrCoeff
        int nCand2 = 1;
        std::partial_sort(Idx2.begin(), Idx2.begin()+nCand2, Idx2.end(),
                          [&](int i, int j){return candidatesFirstLoopCoeff2[i] > candidatesFirstLoopCoeff2[j]; });

        double rotRef, tiltRef;
        for(int i = 0; i < nCand2; i++){

            // reading info of reference image candidate
            mdRef.getRow(rowRef, size_t( candidatesFirstLoop2[ Idx2[i] ] ) );
            rowRef.getValue(MDL_ANGLE_ROT, rotRef);
            rowRef.getValue(MDL_ANGLE_TILT, tiltRef);

            //save metadata of images with angles
            idxOut = mdOut.addObject();
            mdOut.setValue(MDL_IMAGE,       fnImgExp, idxOut);
            mdOut.setValue(MDL_ENABLED,     1, idxOut);
            mdOut.setValue(MDL_IDX,         size_t(candidatesFirstLoop2[ Idx2[i] ]), idxOut);
            mdOut.setValue(MDL_MAXCC,       candidatesFirstLoopCoeff2[Idx2[i]], idxOut);
            mdOut.setValue(MDL_ANGLE_ROT,   rotRef, idxOut);
            mdOut.setValue(MDL_ANGLE_TILT,  tiltRef, idxOut);
            mdOut.setValue(MDL_ANGLE_PSI,    bestRot2[Idx2[i]], idxOut);
            mdOut.setValue(MDL_SHIFT_X,     -1. * bestTx2[Idx2[i]], idxOut);
            mdOut.setValue(MDL_SHIFT_Y,     -1. * bestTy2[Idx2[i]], idxOut);
        }
    }

    // write output metaData file
    mdOut.write(fnDir.getString()+String("outfileMD.xmd"));
    transformerImage.cleanup();
    transformerPolarImage.cleanup();

}

/* print in console some values of double MultidimArray */
void ProgAngularAssignmentMag::printSomeValues(MultidimArray<double> &MDa){
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            std::cout << "val: " << DIRECT_A2D_ELEM(MDa,i,j) << std::endl;
}

void ProgAngularAssignmentMag::produceSideinfo()
{
    mdIn.read(fnIn);
    mdRef.read(fnRef);
}

/* Pearson Coeff */
void ProgAngularAssignmentMag::pearsonCorr(MultidimArray<double> &X, MultidimArray<double> &Y, double &coeff){

    MultidimArray<double>   X2(Ydim,Xdim);
    MultidimArray<double>   Y2(Ydim,Xdim);

    _applyCircularMask(X,X2);
    _applyCircularMask(Y,Y2);
    // covariance
    double X_m, Y_m, X_std, Y_std;
    arithmetic_mean_and_stddev(X2, X_m, X_std);
    arithmetic_mean_and_stddev(Y2, Y_m, Y_std);

    double prod_mean = mean_of_products(X2, Y2);
    double covariace = prod_mean - (X_m * Y_m);

    coeff = covariace / (X_std * Y_std);
}

void ProgAngularAssignmentMag::_applyCircularMask(const MultidimArray<double> &in, MultidimArray<double> &out){

    size_t Cf = (size_t)(Ydim/2.0 + 0.5);
    size_t Cc = (size_t)(Xdim/2.0 + 0.5);
    int pixReduc = 1;
    double rad2 = (Cf - pixReduc) * (Cf - pixReduc);
    double val = 0;

    out.initZeros(Ydim,Xdim);

    for(size_t f = 0; f < Ydim; f++){
        for(size_t c = 0; c < Xdim; c++){
            val = (f-Cf)*(f-Cf) + (c-Cc)*(c-Cc);
            if (val < rad2)
                DIRECT_A2D_ELEM(out, f, c) = DIRECT_A2D_ELEM(in,f,c);
        }
    }

}

/* Arithmetic mean and stdDev for Pearson Coeff */
void ProgAngularAssignmentMag::arithmetic_mean_and_stddev( MultidimArray<double> &data, double &avg, double &stddev ){
    data.computeAvgStdev(avg, stddev);
}

/* Mean of products for Pearson Coeff */
double ProgAngularAssignmentMag::mean_of_products(MultidimArray<double> &data1, MultidimArray<double> &data2){
    double total = 0;
    for (int f = 0; f < Ydim; f++){
        for (int c = 0; c < Xdim; c++){
            total += DIRECT_A2D_ELEM(data1,f,c) * DIRECT_A2D_ELEM(data2,f,c);
        }
    }
    return total/(Xdim*Ydim);
}

/* writing out some data to file with an specified size*/
void ProgAngularAssignmentMag::_writeTestFile(MultidimArray<double> &data, const char* fileName,
                                              size_t nFil, size_t nCol){
    std::ofstream outFile(fileName);
    for (int f = 0; f < nFil; f++){
        for (int c = 0; c < nCol; c++){
            outFile <<  DIRECT_A2D_ELEM(data,f,c) << "\t";
        }
        outFile << "\n";
    }
    outFile.close();
}

/* writing out some data to file Ydim x Xdim size*/
void ProgAngularAssignmentMag::_writeTestFile(MultidimArray<double> &data, const char* fileName){
    std::ofstream outFile(fileName);
    for (int f = 0; f < Ydim; f++){
        for (int c = 0; c < Xdim; c++){
            outFile <<  DIRECT_A2D_ELEM(data,f,c) << "\t";
        }
        outFile << "\n";
    }
    outFile.close();
}

/* get COMPLETE fourier spectrum of Images. It should be changed for half */
void ProgAngularAssignmentMag::_applyFourierImage(MultidimArray<double> &data,
                                                  MultidimArray< std::complex<double> > &FourierData){
    transformerImage.completeFourierTransform(data, FourierData);
}

/* get COMPLETE fourier spectrum of polarRepresentation of Magnitude. It should be changed for half */
void ProgAngularAssignmentMag::_applyFourierImage(MultidimArray<double> &data,
                                                  MultidimArray< std::complex<double> > &FourierData, const size_t &ang){
    transformerPolarImage.completeFourierTransform(data, FourierData);
}

/*first try in using only one half of Fourier space*/
void ProgAngularAssignmentMag::_applyFourierImage2(MultidimArray<double> &data,
                                                   MultidimArray< std::complex<double> > &FourierData){
    transformerImage.FourierTransform(data,FourierData,true);
}

/* first try one half of fourier spectrum of polarRepresentation of Magnitude*/
void ProgAngularAssignmentMag::_applyFourierImage2(MultidimArray<double> &data,
                                                  MultidimArray< std::complex<double> > &FourierData, const size_t &ang){
    transformerPolarImage.FourierTransform(data,FourierData,true); // false --> true para generar copia
}


/* get magnitude of fourier spectrum */
void ProgAngularAssignmentMag::_getComplexMagnitude( MultidimArray< std::complex<double> > &FourierData,
                                                     MultidimArray<double> &FourierMag){
    FFT_magnitude(FourierData,FourierMag);
}

/* cartImg contains cartessian  grid representation of image,
*  rad and ang are the number of radius and angular elements*/
MultidimArray<double> ProgAngularAssignmentMag::imToPolar(MultidimArray<double> &cartIm,
                                                          const size_t &startBand,
                                                          const size_t &finalBand,
                                                          const size_t &n_bands,
                                                          const size_t &rad, const size_t &ang){
    MultidimArray<double> polarImg(n_bands, ang);
    float pi = 3.141592653;
    // coordinates of center
    double cy = (Ydim+1)/2.0;
    double cx = (Xdim+1)/2.0;
    // scale factors
    double sfy = (Ydim-1)/2.0;
    double sfx = (Xdim-1)/2.0;

    double delR = (double)(1.0 / (rad-1));
    double delT = 2.0 * pi / ang;

    // loop through rad and ang coordinates
    double r, t, x_coord, y_coord;
    for(size_t ri = startBand; ri < finalBand; ri++){
        for(size_t ti = 0; ti < ang; ti++ ){
            r = ri * delR;
            t = ti * delT;
            x_coord = ( r * cos(t) ) * sfx + cx;
            y_coord = ( r * sin(t) ) * sfy + cy;
            // set value of polar img
            DIRECT_A2D_ELEM(polarImg,ri-startBand,ti) = interpolate(cartIm,x_coord,y_coord);
        }
    }

    return polarImg;
}

/* cartImg contains cartessian  grid representation of image,
*  rad and ang are the number of radius and angular elements
*  this function was built for half representation of Fourier spectrum*/
MultidimArray<double> ProgAngularAssignmentMag::imToPolar2(MultidimArray<double> &cartIm,
                                                          const size_t &rad, const size_t &ang){
    MultidimArray<double> polarImg(rad, ang);
    float pi = 3.141592653;
    // coordinates of center
    double cy = 0.5; //(Ydim+1)/2.0;
    double cx = (Xdim+1)/2.0;
    // scale factors
    double sfy = (Ydim-1)/2.0;
    double sfx = (Xdim-1)/2.0;

    double delR = (double)(1.0 / (rad-1));
    double delT = pi / ang;

    // loop through rad and ang coordinates
    double r, t, x_coord, y_coord;
    for(size_t ri = 0; ri < rad; ri++){
        for(size_t ti = 0; ti < ang; ti++ ){
            r = ri * delR;
            t = ti * delT;
            x_coord = ( r * cos(t) ) * sfx + cx;
            y_coord = ( r * sin(t) ) * sfy + cy;

            // set value of polar img
            DIRECT_A2D_ELEM(polarImg,ri,ti) = interpolate(cartIm,x_coord,y_coord);
        }
    }


    return polarImg;
}

/* bilinear interpolation */
double ProgAngularAssignmentMag::interpolate(MultidimArray<double> &cartIm,
                                             double &x_coord, double &y_coord){
    double val;
    size_t xf = floor(x_coord);
    size_t xc = ceil(x_coord);
    size_t yf = floor(y_coord);
    size_t yc = ceil(y_coord);

    if ( (xf == xc) && ( yf == yc )){
        val = dAij(cartIm, xc, yc);
    }
    else if (xf == xc){ // linear
        val = dAij(cartIm, xf, yf) + (y_coord - yf) * ( dAij(cartIm, xf, yc) - dAij(cartIm, xf, yf) );
    }
    else if(yf == yc){ // linear
        val = dAij(cartIm, xf, yf) + (x_coord - xf) * ( dAij(cartIm, xc, yf) - dAij(cartIm, xf, yf) );
    }
    else{ // bilinear
        val = ((double)(( dAij(cartIm,xf,yf)*(yc-y_coord) + dAij(cartIm,xf,yc)*(y_coord-yf) ) * (xc - x_coord)) +
               (double)(( dAij(cartIm,xc,yf)*(yc-y_coord) + dAij(cartIm,xc,yc)*(y_coord-yf) ) * (x_coord - xf))
               )  / (double)( (xc - xf)*(yc - yf) );
    }

    return val;

}

/* its an experiment for implement fftshift*/
void ProgAngularAssignmentMag::completeFourierShift(MultidimArray<double> &in, MultidimArray<double> &out){

    // correct output size
    out.resizeNoCopy(in);

    size_t Cf = (size_t)(YSIZE(in)/2.0 + 0.5);      //(Ydim/2.0 + 0.5);
    size_t Cc = (size_t)(XSIZE(in)/2.0 + 0.5);      //(Xdim/2.0 + 0.5);

    size_t ff, cc;
    for(size_t f = 0; f < YSIZE(in); f++){
        ff = (f + Cf) % YSIZE(in);
        for(size_t c = 0; c < XSIZE(in); c++){
            cc = (c + Cc) % XSIZE(in);
            DIRECT_A2D_ELEM(out, ff, cc) = DIRECT_A2D_ELEM(in,f,c);
        }
    }
}

/* its an experiment for implement fftshift*/
void ProgAngularAssignmentMag::halfFourierShift(MultidimArray<double> &in, MultidimArray<double> &out){
    size_t Cf = (size_t)(Ydim/2.0 + 0.5);
    out.resizeNoCopy(in);

    size_t ff, cc;
    for(size_t f = 0; f < Ydim; f++){
        ff = (f + Cf) % Ydim;
        for(size_t c = 0; c < Cf; c++){
            cc = c;
            DIRECT_A2D_ELEM(out, ff, cc) = DIRECT_A2D_ELEM(in,f,c);
        }
    }

}



/* experiment for GCC matrix product F1 .* conj(F2)
*
*/
void ProgAngularAssignmentMag::ccMatrix(MultidimArray< std::complex<double>> &F1,
                                        MultidimArray< std::complex<double>> &F2,
                                        MultidimArray<double> &result){


    result.resizeNoCopy(YSIZE(F1),2*(XSIZE(F1)-1));

    CorrelationAux aux;
    aux.transformer1.setReal(result);
    aux.transformer1.setFourier(F1);
    // Multiply FFT1 * FFT2'
    double a, b, c, d; // a+bi, c+di
    double *ptrFFT2=(double*)MULTIDIM_ARRAY(F2);
    double *ptrFFT1=(double*)MULTIDIM_ARRAY(aux.transformer1.fFourier);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
    {
        a=*ptrFFT1;
        b=*(ptrFFT1+1);
        c=(*ptrFFT2++);
        d=(*ptrFFT2++)*(-1);
        *ptrFFT1++ = a*c-b*d;
        *ptrFFT1++ = b*c+a*d;
    }

    aux.transformer1.inverseFourierTransform();
    CenterFFT(result, true);
    result.setXmippOrigin();
}


/* select n_bands of polar representation of magnitude spectrum */
void ProgAngularAssignmentMag::selectBands(MultidimArray<double> &in, MultidimArray<double> &out,
                                           const size_t &n_bands, const size_t &startBand, const size_t &n_ang ){

    int colStop = XSIZE(out);
    int rowStop = YSIZE(out);
    int i, j;
    // 0:179 and duplicate data
    for (i = 0; i < rowStop; i++){
        for (j = 0; j < colStop; j++){
            dAij(out,i,j) = dAij(in, startBand+i, j);
        }
    }

}

/* gets maximum value for each column*/
void ProgAngularAssignmentMag::maxByColumn(MultidimArray<double> &in,
                                           MultidimArray<double> &out,
                                           const size_t &nFil, const size_t &nCol){

    out.resizeNoCopy(1,XSIZE(in));
    int f, c;
    double maxVal, val2;
    for(c = 0; c < XSIZE(in); c++){
        maxVal = dAij(in, 0, c);
        for(f = 1; f < YSIZE(in); f++){
            val2 = dAij(in, f, c);
            if (val2 > maxVal)
                maxVal = val2;
        }
        dAi(out,c) = maxVal;
    }
}

/* gets maximum value for each row */
void ProgAngularAssignmentMag::maxByRow(MultidimArray<double> &in,
                                        MultidimArray<double> &out,
                                        const size_t &nFil, const size_t &nCol){
    out.resizeNoCopy(1,YSIZE(in));
    int f, c;
    double maxVal, val2;
    for(f = 0; f < nFil; f++){
        maxVal = dAij(in, f, 0);
        for(c = 1; c < nCol; c++){
            val2 = dAij(in, f, c);
            if (val2 > maxVal)
                maxVal = val2;
        }
        dAi(out,f) = maxVal;
    }
}

/*quadratic interpolation for location of peak in crossCorr vector*/
double quadInterp(const int Idx, MultidimArray<double> &in){
    double InterpIdx = Idx - ( ( dAi(in,Idx+1) - dAi(in,Idx-1) ) / ( dAi(in,Idx+1) + dAi(in,Idx-1) - 2*dAi(in, Idx) ) )/2.;
    return InterpIdx;
}

/* Only for 180 angles */
/* approach which selects only two locations of maximum peaks in ccvRot */
void ProgAngularAssignmentMag::rotCandidates3(MultidimArray<double> &in,
                                              std::vector<double> &cand,
                                              const size_t &size, int *nPeaksFound){
    const int maxNumPeaks = 20;
    double max1 = -10.;
    int idx1 = 0;
    double max2 = -10.;
    int idx2 = 0;
    int i;
    int cont = 0;
    *(nPeaksFound) = cont;

    for(i = 89/*1*/; i < 271/*size-1*/; i++){ // only look for in range -90:90
        // current value is a peak value?
        if ( (dAi(in,i) > dAi(in,i-1)) && (dAi(in,i) > dAi(in,i+1)) ){
            cont++;
            if( dAi(in,i) > max1){
                max2 = max1;
                idx2 = idx1;
                max1 = dAi(in,i);
                idx1 = i;
            }
            else if( dAi(in,i) > max2 && dAi(in,i) != max1 ){
                max2 = dAi(in,i);
                idx2 = i;
            }
        }
    }

    if( cont > maxNumPeaks){
        printf("reaches max number of peaks!\n");

    }

    int maxAccepted = 2;

    maxAccepted = ( cont < maxAccepted) ? cont : maxAccepted;

    if(cont){
        std::vector<int> temp(2,0);
        temp[0] = idx1;
        temp[1] = idx2;
        int tam = 2*maxAccepted;
        *(nPeaksFound) = tam;
        cand.reserve(tam);
        double interpIdx; // quadratic interpolated location of peak
        double interpDiff;
        for(i = 0; i < maxAccepted; i++){
            interpIdx = quadInterp(temp[i], in);
            interpDiff = temp[i] - interpIdx;
            cand[i] =  double( size - 1 )/2. - interpIdx;
            cand[i+maxAccepted] =(cand[i]>0) ? cand[i] + 180 : cand[i] - 180 ;
        }

    }
    else{
        printf("no peaks found!\n");
    }
}

/*method in xmipp_fftw but here I make a copy of FFT1*/
void ProgAngularAssignmentMag::fast_correlation_vector2(MultidimArray<std::complex<double> > FFT1,
                                                        const MultidimArray<std::complex<double> > FFT2,
                                                        MultidimArray<double> &R,
                                                        FourierTransformer &transformer){


}

/* approach which selects only two locations of maximum peaks in ccvRot */
void ProgAngularAssignmentMag::rotCandidates2(MultidimArray<double> &in,
                                              std::vector<double> &cand,
                                              const size_t &size, int *nPeaksFound){
    const int maxNumPeaks = 20;
    double max1 = -10.;
    int idx1 = 0;
    double max2 = -10.;
    int idx2 = 0;
    int i;
    int cont = 0;
    *(nPeaksFound) = cont;
    for(i = 89/*1*/; i < 271/*size-1*/; i++){
        // current value is a peak value?
        if ( (dAi(in,i) > dAi(in,i-1)) && (dAi(in,i) > dAi(in,i+1)) ){
            cont++;
            if( dAi(in,i) > max1){
                max2 = max1;
                idx2 = idx1;
                max1 = dAi(in,i);
                idx1 = i;
            }
            else if( dAi(in,i) > max2 && dAi(in,i) != max1 ){
                max2 = dAi(in,i);
                idx2 = i;
            }
        }
    }

    if( cont > maxNumPeaks){
        printf("reaches max number of peaks!\n");
    }

    int maxAccepted = 2;

    maxAccepted = ( cont < maxAccepted) ? cont : maxAccepted;

    if(cont){
        std::vector<int> temp(2,0);
        temp[0] = idx1;
        temp[1] = idx2;
        int tam = 2*maxAccepted;
        *(nPeaksFound) = tam;
        cand.reserve(tam);
        for(i = 0; i < maxAccepted; i++){
            cand[i] = dAi(axRot,temp[i]);
            cand[i+maxAccepted] =(cand[i]>0) ? cand[i] + 180 : cand[i] - 180 ;
        }

    }
    else{
        printf("no peaks found!\n");
    }
}

/* candidates to best rotation*/
void ProgAngularAssignmentMag::rotCandidates(MultidimArray<double> &in,
                                             std::vector<double> &cand,
                                             const size_t &size, int *nPeaksFound){
    const int maxNumPeaks = 30;
    int maxAccepted = 4;
    int *peakPos = (int*) calloc(maxNumPeaks,sizeof(int));
    int cont = 0;
    *(nPeaksFound) = cont;
    int i;
    for(i = 89/*1*/; i < 271/*size-1*/; i++){

        if ( (dAi(in,i) > dAi(in,i-1)) && (dAi(in,i) > dAi(in,i+1)) ){
            peakPos[cont] = i;
            cont++;
            *(nPeaksFound) = cont;
        }
    }

    maxAccepted = ( *(nPeaksFound) < maxAccepted) ? *(nPeaksFound) : maxAccepted;

    if( *(nPeaksFound) > maxNumPeaks)
        printf("reaches max number of peaks!\n");

    if(cont){
        std::vector<int> temp(*(nPeaksFound),0);
        for(i = 0; i < *(nPeaksFound); i++){
            temp[i] = peakPos[i];
        }
        // delete peakPos
        free(peakPos);

        // sorting first in case there are more than maxAccepted peaks
        std::sort(temp.begin(), temp.end(), [&](int i, int j){return dAi(in,i) > dAi(in,j); } );

        int tam = 2*maxAccepted; //
        *(nPeaksFound) = tam;
        cand.reserve(tam);
        for(i = 0; i < maxAccepted; i++){
            cand[i] = dAi(axRot,temp[i]);
            cand[i+maxAccepted] =(cand[i]>0) ? cand[i] + 180 : cand[i] - 180 ;
        }
    }
    else{
        printf("no peaks found!\n");
        // delete peakPos
        free(peakPos);
    }

}

/* instace of "delay axes" for assign rotation and traslation candidates*/
void ProgAngularAssignmentMag::_delayAxes(const size_t &Ydim, const size_t &Xdim, const size_t &n_ang){
    axRot.resize(1,1,1,n_ang);
    axTx.resize(1,1,1,Xdim);
    axTy.resize(1,1,1,Ydim);

    double M = double(n_ang - 1)/2.;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(axRot){
        dAi(axRot,i) = ceil(M - i);
    }
    M = double(Xdim - 1)/2.0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(axTx){
        dAi(axTx,i) = ceil(M - i);
    }
    M = double(Ydim - 1)/2.0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(axTy){
        dAi(axTy,i) = ceil(M - i);
    }
}

/* selection of best candidate to rotation and its corresponding shift
 * called at first loop in "coarse" searching
 * shitfs are computed as maximum of CrossCorr vector
 * vector<double> cand contains candidates to relative rotation between images
*/
void ProgAngularAssignmentMag::bestCand(/*inputs*/
                                        MultidimArray<double> &MDaIn,
                                        MultidimArray< std::complex<double> > &MDaInF,
                                        MultidimArray<double> &MDaRef,
                                        std::vector<double> &cand,
                                        int &peaksFound,
                                        /*outputs*/
                                        double *bestCandRot,
                                        double *shift_x,
                                        double *shift_y,
                                        double *bestCoeff){
    *(bestCandRot) = 0;
    *(shift_x) = 0.;
    *(shift_y) = 0.;
    *(bestCoeff) = -10.0;
    double rotVar = 0.0;
    double tempCoeff;
    double tx, ty;
    MultidimArray<double> MDaRefRot;
    MultidimArray<double> MDaRefRotShift;
    MultidimArray<double> ccMatrixShift;
    MultidimArray<double> ccVectorTx;
    MultidimArray<double> ccVectorTy;
    MultidimArray< std::complex<double> > MDaRefRotF;

    MDaRefRot.setXmippOrigin();
    for(int i = 0; i < peaksFound; i++){
        rotVar = -1. * cand[i]; //
        _applyRotation(MDaRef,rotVar,MDaRefRot);

        _applyFourierImage2(MDaRefRot,MDaRefRotF);

        ccMatrix(MDaInF, MDaRefRotF, ccMatrixShift);// cross-correlation matrix
        maxByColumn(ccMatrixShift, ccVectorTx, YSIZE(ccMatrixShift), XSIZE(ccMatrixShift)); // ccvMatrix to ccVector
        getShift(axTx, ccVectorTx,tx,XSIZE(ccMatrixShift));
        tx = -1. * tx;
        maxByRow(ccMatrixShift, ccVectorTy, YSIZE(ccMatrixShift), XSIZE(ccMatrixShift)); // ccvMatrix to ccVector
        getShift(axTy, ccVectorTy,ty,YSIZE(ccMatrixShift));
        ty = -1. * ty;

        if ( std::abs(tx)>10 || std::abs(ty)>10 ) // 10 es elegido pero debo poner criterio automático, por ejemplo Xdim*0.1
            continue;

        _applyShift(MDaRefRot, tx, ty, MDaRefRotShift);

        // Pearson coeff
        pearsonCorr(MDaIn, MDaRefRotShift, tempCoeff);

        if ( tempCoeff > *(bestCoeff) ){
            *(bestCoeff) = tempCoeff;
            *(shift_x) = tx;
            *(shift_y) = ty;
            *(bestCandRot) = rotVar;
        }

    }

}

/* apply rotation */
void ProgAngularAssignmentMag::_applyRotation(MultidimArray<double> &MDaRef, double &rot,
                                              MultidimArray<double> &MDaRefRot){
    // Transform matrix
    Matrix2D<double> A(3,3);
    A.initIdentity();
    double ang, cosine, sine;
    ang = DEG2RAD(rot);
    cosine = cos(ang);
    sine = sin(ang);

    // rotation
    MAT_ELEM(A,0, 0) = cosine;
    MAT_ELEM(A,0, 1) = sine;
    MAT_ELEM(A,1, 0) = -sine;
    MAT_ELEM(A,1, 1) = cosine;

    // Shift
    MAT_ELEM(A,0, 2) = 0.;
    MAT_ELEM(A,1, 2) = 0.;

    applyGeometry(LINEAR, MDaRefRot, MDaRef, A, IS_NOT_INV, DONT_WRAP);

}

/* apply traslation */
void ProgAngularAssignmentMag::_applyShift(MultidimArray<double> &MDaRef,
                                           double &tx, double &ty,
                                           MultidimArray<double> &MDaRefShift){
    // Transform matrix
    Matrix2D<double> A(3,3);
    A.initIdentity();

    // Shift
    MAT_ELEM(A,0, 2) = tx;
    MAT_ELEM(A,1, 2) = ty;

    applyGeometry(LINEAR, MDaRefShift, MDaRef, A, IS_NOT_INV, DONT_WRAP);
}

/* finds shift as maximum of ccVector */
void ProgAngularAssignmentMag::getShift(MultidimArray<double> &axis,
                                        MultidimArray<double> &ccVector, double &shift, const size_t &size){
    double maxVal = -10.;
    int idx;
    int i;
    for(i = 0; i < size; i++){
        if(ccVector[i] > maxVal){
            maxVal = ccVector[i];
            idx = i;
        }
    }

    //    shift = dAi(axis, idx);

    // interpolate value
    double interpIdx, interpDiff;
    interpIdx = quadInterp(idx, ccVector);
    interpDiff = idx - interpIdx;
    shift = double( size - 1 )/2. - interpIdx;

}


/* Structural similarity SSIM index Coeff */
void ProgAngularAssignmentMag::ssimIndex(MultidimArray<double> &X, MultidimArray<double> &Y, double &coeff){

    // covariance
    double X_m, Y_m, X_std, Y_std;
    double c1, c2, L;
    arithmetic_mean_and_stddev(X, X_m, X_std);
    arithmetic_mean_and_stddev(Y, Y_m, Y_std);

    double prod_mean = mean_of_products(X, Y);
    double covariace = prod_mean - (X_m * Y_m);

    L = 1;
    c1 = (0.01*L) * (0.01*L);
    c2 = (0.03*L) * (0.03*L); // estabilidad en división


    coeff = ( (2*X_m*Y_m + c1)*(2*covariace+c2) )/( (X_m*X_m + Y_m*Y_m + c1)*(X_std*X_std + Y_std*Y_std + c2) );
}

/* selection of best candidate to rotation and its corresponding shift
 * called at second loop in a little bit more strict searching
 * shitfs are computed as maximum of CrossCorr vector +0.5 / -0.5
 * vector<double> cand contains candidates to relative rotation between images with larger CrossCorr-coeff after first loop
*/
void ProgAngularAssignmentMag::bestCand2(/*inputs*/
                                         MultidimArray<double> &MDaIn,
                                         MultidimArray< std::complex<double> > &MDaInF,
                                         MultidimArray<double> &MDaRef,
                                         std::vector<double> &cand,
                                         int &peaksFound,
                                         /*outputs*/
                                         double *bestCandRot,
                                         double *shift_x,
                                         double *shift_y,
                                         double *bestCoeff){
    *(bestCandRot) = 0;
    *(shift_x) = 0.;
    *(shift_y) = 0.;
    *(bestCoeff) = -10.0;
    double rotVar = 0.0;
    double tempCoeff;
    double tx, ty;
    std::vector<double> vTx, vTy;
    MultidimArray<double> MDaRefRot;
    MultidimArray<double> MDaRefRotShift;
    MultidimArray<double> ccMatrixShift;
    MultidimArray<double> ccVectorTx( (const size_t) 1,XSIZE(MDaInF));
    MultidimArray<double> ccVectorTy( (const size_t) 1,YSIZE(MDaInF));
    MultidimArray< std::complex<double> > MDaRefRotF;

    MDaRefRot.setXmippOrigin();
    for(int i = 0; i < peaksFound; i++){
        rotVar = -1. * cand[i];
        _applyRotation(MDaRef,rotVar,MDaRefRot);

        _applyFourierImage2(MDaRefRot,MDaRefRotF); // fourier --> F2_r

        ccMatrix(MDaInF, MDaRefRotF, ccMatrixShift);// cross-correlation matrix
        maxByColumn(ccMatrixShift, ccVectorTx, YSIZE(ccMatrixShift), XSIZE(ccMatrixShift)); // ccvMatrix to ccVector
        getShift(axTx, ccVectorTx,tx,XSIZE(ccMatrixShift));
        tx = -1. * tx;
        maxByRow(ccMatrixShift, ccVectorTy, YSIZE(ccMatrixShift), XSIZE(ccMatrixShift)); // ccvMatrix to ccVector
        getShift(axTy, ccVectorTy,ty,YSIZE(ccMatrixShift));
        ty = -1. * ty;

        if ( std::abs(tx)>10 || std::abs(ty)>10 ) // 10 es elegido pero debo poner criterio automático
            continue;

        //*********** when strict, after first loop ***************
        // posible shifts
        vTx.push_back(tx);
        vTx.push_back(tx+0.5);
        vTx.push_back(tx-0.5);
        vTy.push_back(ty);
        vTy.push_back(ty+0.5);
        vTy.push_back(ty-0.5);

        for(int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                // translate rotated version of MDaRef
                _applyShift(MDaRefRot, vTx[j], vTy[k], MDaRefRotShift);
                // Pearson coeff

                pearsonCorr(MDaIn, MDaRefRotShift, tempCoeff);
                //        std::cout << "myCorr(f1,f2_rt): " << tempCoef << std::endl;
                if ( tempCoeff > *(bestCoeff) ){
                    *(bestCoeff) = tempCoeff;
                    *(shift_x) = vTx[j];
                    *(shift_y) = vTy[k];
                    *(bestCandRot) = rotVar;
                }
            }
        }


    }

}

/* apply rotation */
void ProgAngularAssignmentMag::_applyRotationAndShift(MultidimArray<double> &MDaRef, double &rot, double &tx, double &ty,
                                                      MultidimArray<double> &MDaRefRot){
    // Transform matrix
    Matrix2D<double> A(3,3);
    A.initIdentity();
    double ang, cosine, sine;
    ang = DEG2RAD(rot);
    cosine = cos(ang);
    sine = sin(ang);

    // rotation
    MAT_ELEM(A,0, 0) = cosine;
    MAT_ELEM(A,0, 1) = sine;
    MAT_ELEM(A,1, 0) = -sine;
    MAT_ELEM(A,1, 1) = cosine;

    // Shift
    MAT_ELEM(A,0, 2) = tx;
    MAT_ELEM(A,1, 2) = ty;

    applyGeometry(LINEAR, MDaRefRot, MDaRef, A, IS_NOT_INV, DONT_WRAP);

}



