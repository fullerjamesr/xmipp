/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

#ifndef __ANGULAR_ASSIGNMENT_MAG_H
#define __ANGULAR_ASSIGNMENT_MAG_H

#include <../../xmippCore/core/xmipp_program.h>
#include <../../xmippCore/core/xmipp_fftw.h>
#include <../../xmippCore/core/metadata_extension.h>
#include <../../xmippCore/core/multidim_array.h>
#include <../../xmipp/libraries/data/mask.h>

#include <vector>
#include <fstream> // borrar luego
#include <ctime>


/**@defgroup AngularAssignmentMag ***
   @ingroup ReconsLibrary */
//@{

/** Angular Assignment mag parameters. */
class ProgAngularAssignmentMag: public XmippProgram
{
public:
    /** Filenames */
    FileName fnIn, fnDir, fnSym, fnRef;
public: // Internal members
    // Metadata with input images and input volumes
    MetaData mdIn, mdRef, mdOut;

    // vector of reference images
    std::vector< MultidimArray<double> > vecMDaRef;

    // vector of Fourier of reference images
    std::vector< MultidimArray< std::complex<double> > > vecMDaRefF;

    // vector of Fourier of polar representation of magnitude spectrum of reference images
    std::vector< MultidimArray< std::complex<double> > > vecMDaRefFMs_polarF;

    // Size of the images
    size_t Xdim, Ydim;

    // Transformer
    FourierTransformer transformerImage, transformerPolarImage;


    // "delay axes"
    MultidimArray<double> axRot;
    MultidimArray<double> axTx;
    MultidimArray<double> axTy;

public:
    /// Read arguments from command line
    void defineParams();
    void readParams();

    /** Show. */
    void show();

    /** Run. */
    void run();

    /// Produce side info: fill arrays with relevant transformation matrices
    void produceSideinfo();

private:
    void printSomeValues(MultidimArray<double> & MDa);
    void pearsonCorr(MultidimArray<double> &X, MultidimArray<double> &Y, double &coeff);
    void arithmetic_mean_and_stddev(MultidimArray<double> &data, double &avg, double &stddev);
    double mean_of_products(MultidimArray<double> &data1, MultidimArray<double> &data2);
    void _writeTestFile(MultidimArray<double> &data, const char *fileName);
    void _writeTestFile(MultidimArray<double> &data, const char *fileName, size_t nFil, size_t nCol);
    void _applyFourierImage(MultidimArray<double> &data, MultidimArray<std::complex<double> > &FourierData);
    void _applyFourierImage(MultidimArray<double> &data, MultidimArray<std::complex<double> > &FourierData, const size_t &ang);
    void _getComplexMagnitude(MultidimArray<std::complex<double> > &FourierData, MultidimArray<double> &FourierMag);
    MultidimArray<double> imToPolar(MultidimArray<double> &cartIm, const size_t &rad, const size_t &ang);
    double interpolate(MultidimArray<double> &cartIm, double &x_coord, double &y_coord);
    void completeFourierShift(MultidimArray<double> &in, MultidimArray<double> &out);
    void ccMatrix(MultidimArray<std::complex<double> > &F1, MultidimArray<std::complex<double> > &F2, MultidimArray<double> &result);
    void selectBands(MultidimArray<double> &in, MultidimArray<double> &out, const size_t &n_bands, const size_t &startBand, const size_t &n_ang);
    void maxByColumn(MultidimArray<double> &in, MultidimArray<double> &out, const size_t &nFil, const size_t &nCol);
    void rotCandidates(MultidimArray<double> &in, std::vector<double>& cand, const size_t &size, int *nPeaksFound);
    void _delayAxes(const size_t &Ydim, const size_t &Xdim, const size_t &n_ang);
    void bestCand(MultidimArray<double> &MDaIn, MultidimArray<std::complex<double> > &MDaInF, MultidimArray<double> &MDaRef, std::vector<double> &cand, int &peaksFound, double *bestCandRot, double *shift_x, double *shift_y, double *bestCoeff);
    void _applyRotation(MultidimArray<double> &MDaRef, double &rot, MultidimArray<double> &MDaRefRot);
    void maxByRow(MultidimArray<double> &in, MultidimArray<double> &out, const size_t &nFil, const size_t &nCol);
    void getShift(MultidimArray<double> &axis, MultidimArray<double> &ccVector, double &shift, const size_t &size);
    void _applyShift(MultidimArray<double> &MDaRef, double &tx, double &ty, MultidimArray<double> &MDaRefShift);
    void ssimIndex(MultidimArray<double> &X, MultidimArray<double> &Y, double &coeff);
    void bestCand2(MultidimArray<double> &MDaIn, MultidimArray<std::complex<double> > &MDaInF, MultidimArray<double> &MDaRef, std::vector<double> &cand, int &peaksFound, double *bestCandRot, double *shift_x, double *shift_y, double *bestCoeff);
    void _applyRotationAndShift(MultidimArray<double> &MDaRef, double &rot, double &tx, double &ty, MultidimArray<double> &MDaRefRot);
    void rotCandidates2(MultidimArray<double> &in, std::vector<double> &cand, const size_t &size, int *nPeaksFound);
    void _applyFourierImage2(MultidimArray<double> &data, MultidimArray<std::complex<double> > &FourierData);
    void _applyFourierImage2(MultidimArray<double> &data, MultidimArray<std::complex<double> > &FourierData, const size_t &ang);
    void halfFourierShift(MultidimArray<double> &in, MultidimArray<double> &out);
    MultidimArray<double> imToPolar2(MultidimArray<double> &cartIm, const size_t &rad, const size_t &ang);
    void rotCandidates3(MultidimArray<double> &in, std::vector<double> &cand, const size_t &size, int *nPeaksFound);
    void fast_correlation_vector2(MultidimArray<std::complex<double> > FFT1, const MultidimArray<std::complex<double> > FFT2, MultidimArray< double >& R, FourierTransformer &transformer);

    void _applyCircularMask(const MultidimArray<double> &in, MultidimArray<double> &out);
};
//@}



//            // test  write MdaRef
//            _writeTestFile(MDaRef,"/home/jeison/Escritorio/t_ref.txt");

//            // test fftw and magnitude
//            MultidimArray< std::complex<double> > MDaRef_fourier(Ydim, Xdim);
//            _applyFourier(MDaRef, MDaRef_fourier);
//            MultidimArray<double> MDaRef_fourier_mag(Ydim, Xdim);
//            _getComplexMagnitude(MDaRef_fourier,MDaRef_fourier_mag);
//            _writeTestFile(MDaRef_fourier_mag,"/home/jeison/Escritorio/t_mag.txt");

            // test  write MDaRef after fourier
            // supongo que no tengo problema porque genero todo el espectro
            // genero una copia y se supone que así es más lento
//            _writeTestFile(MDaRef,"/home/jeison/Escritorio/t_refAterFourier.txt");

//            // test magnitude fourier shift
//            MultidimArray< double > MDaRef_fourier_mag_shifted(Ydim, Xdim);
//            completeFourierShift(MDaRef_fourier_mag, MDaRef_fourier_mag_shifted);
//            _writeTestFile(MDaRef_fourier_mag_shifted,"/home/jeison/Escritorio/t_magShifted.txt");

//            // test imTopolar
//            const size_t n_rad = size_t(Xdim/2 + 0.5);
//            const size_t n_ang = size_t(360);
//            MultidimArray<double> MDaRef_polar(n_rad,n_ang);
//            MDaRef_polar = imToPolar(MDaRef, n_rad, n_ang);
//            // test write output
//            _writeTestFile(MDaRef_polar,"/home/jeison/Escritorio/t_polar.txt", n_rad, n_ang);

//            // pearson test
//            double coeff;
//            pearsonCorr(MDaIn, MDaRef,coeff);
//            std::cout << "coeff: " << coeff << std::endl;

//                // test de organización
//                int ints[] = {0,1,2,3,4,5};
//                std::vector<int> idx(ints,ints + sizeof(ints)/sizeof(int));
//                double ints2[] = {1.3,2.5,5.2,8.3,16.6,4.3};
//                std::vector<double> Ordidx(ints2,ints2 + sizeof(ints2)/sizeof(double));
//                std::cout << "before order\n";
//                for(int k = 0; k < 6; k++)
//                    std::cout << "idx("<<k<<")= "<<idx[k]<<"\t OrdIdx("<<k<<")= "<<Ordidx[k]<<std::endl;
//                std::partial_sort(idx.begin(), idx.begin()+3, idx.end(), [&](int i, int j){return ints2[i] > ints2[j]; });
//                std::cout << "after partial order\n";
//                for(int k = 0; k < 6; k++)
//                    std::cout << "idx("<<k<<")= "<<idx[k]<<"\t OrdIdx("<<k<<")= "<<Ordidx[k]<<std::endl;



//        // show first loop results
//        //std::cout << "input image " << countInImg + 1 << ", direction candidates: \n";
//        outfile << "size of IdxVector underwent to sorting: " << Idx.size() << "\n";
//        outfile << "input image " << countInImg + 1 << ", direction candidates: \n";

//        for(int i = 0; i < nCand; i++)
//            outfile   << "dir:  "          << candidatesFirstLoop[ Idx[i] ]
//                      << "\t    coef:  "    << candidatesFirstLoopCoeff[Idx[i]]
//                      << "\t    rot:  "     << bestRot[Idx[i]]
//                      << "\t    tx:  "      << bestTx[Idx[i]]
//                      << "\t    ty:  "      << bestTy[Idx[i]] << "\n";
//        outfile << "\t \t real Rot/psi: "   << psiVal
//                << "\t    realTx: "         << realTx
//                << "\t    realTy: "         << realTy << "\n";
//        outfile << "\n";
// clean candidates
//        std::vector<double>().swap(candidatesFirstLoopCoeff);
//        std::vector<unsigned int>().swap(candidatesFirstLoop);

//            // en el ciclo debo elegir los mejores candidatos (anterior)
//            // antes de incluir partial_sort
//            double thres = tempCoeff - 0.02; //before tempCoeff * (1. - 0.03); // 0.02 -- 0.01
//            if(bestCoeff > thres){
//                tempCoeff = bestCoeff;
//                Idx.push_back(k++);
//                candidatesFirstLoop.push_back(countRefImg+1);
//                candidatesFirstLoopCoeff.push_back(bestCoeff);
//                bestTx.push_back(Tx);
//                bestTy.push_back(Ty);
//                bestRot.push_back(bestCandVar);
//            }

//            // lectura de imagenes referencia dentro del primer ciclo
//            mdRef.getRow(rowRef, size_t(countRefImg+1) /*iterRef->objId*/);
//            rowRef.getValue(MDL_IMAGE, fnImgRef);
//            // processing reference image
//            ImgRef.read(fnImgRef);
//            MDaRef = ImgRef();
//            _applyFourier(MDaRef, MDaRefF);// fourier experimental image (genera copia?)
//            _getComplexMagnitude(MDaRefF, MDaRefFM);// magnitude espectra experimental image
//            completeFourierShift(MDaRefFM, MDaRefFMs);// shift spectrum
//            MDaRefFMs_polar = imToPolar(MDaRefFMs, n_rad, n_ang);// polar representation of magnitude
//            selectBands(MDaRefFMs_polar, MDaRefFMs_polarPart, n_bands, startBand, n_ang); // select bands
//            _applyFourier(MDaRefFMs_polarPart,MDaRefFMs_polarF); // apply fourier

//            // lectura de imágenes referencia en el segundo ciclo
//            mdRef.getRow(rowRef, size_t(candidatesFirstLoop[ Idx[i] ]));
//            rowRef.getValue(MDL_IMAGE, fnImgRef);
//            // En lugar de hacer todo esto podria almacenar los candidatos a rotación del primer loop
//            // processing reference image
//            ImgRef.read(fnImgRef);
//            MDaRef = ImgRef();
//            // aplicar dicha rotación a la imagen referencia y volver a calcular rotación y traslación
//            double rotVal = bestRot[ Idx[i] ];
//            double trasXval = bestTx[ Idx[i] ];
//            double trasYval = bestTy[ Idx[i] ];
//            _applyRotationAndShift(MDaRef, rotVal, trasXval, trasYval, MDaRefTrans);
//            _applyFourier(MDaRefTrans, MDaRefF);// fourier experimental image (genera copia?)
//            _getComplexMagnitude(MDaRefF, MDaRefFM);// magnitude espectra experimental image
//            completeFourierShift(MDaRefFM, MDaRefFMs);// shift spectrum
//            MDaRefFMs_polar = imToPolar(MDaRefFMs, n_rad, n_ang);// polar representation of magnitude
//            selectBands(MDaRefFMs_polar, MDaRefFMs_polarPart, n_bands, startBand, n_ang); // select bands
//            _applyFourier(MDaRefFMs_polarPart,MDaRefFMs_polarF); // apply fourier

#endif
