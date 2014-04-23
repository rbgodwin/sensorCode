#include "aquila/source/WaveFile.h"
#include <iostream>
#include "aquila/transform/FftFactory.h"
#include "aquila/tools/TextPlot.h"
#include "aquila/tools/TextPlot.h"

#undef DEBUG

int main(int argc, char *argv[])
{
  Aquila::TextPlot plot;
  //  int SIZE=512;
  int SIZE=16384;

  double wavRightEnergy=0;
  double wavLeftEnergy=0;
  double deltaT, deltaTotal=0, deltaE=0;
  int i,j,k;
  int frames;
  int maxCorLoc;
  Aquila::SampleType maxCorVal;
  Aquila::SampleType autoCor[SIZE];
  Aquila::SpectrumType correlation;
  Aquila::SampleType frameLeft[SIZE];
  Aquila::SampleType frameRight[SIZE];
  Aquila::SpectrumType spectrumRight,spectrumLeft;    
  correlation.resize(SIZE);

  auto fft = Aquila::FftFactory::getFft(SIZE);

    if (argc < 2)
    {
        std::cout << "Usage: wave_info <FILENAME>" << std::endl;
        return 1;
    }
    Aquila::WaveFile wav(argv[1],Aquila::LEFT);
    Aquila::WaveFile wavR(argv[1],Aquila::RIGHT);

    std::cout << "Filenames: "           << wav.getFilename();
    std::cout << "\nLength: "           << wav.getAudioLength()     << " ms";
    std::cout << "\nSample frequency: " << wav.getSampleFrequency() << " Hz";
    std::cout << "\nChannels: "         << wav.getChannelsNum();
    std::cout << "\nByte rate: "        << wav.getBytesPerSec()/1024 << " kB/s";
    std::cout << "\nBits per sample: "  << wav.getBitsPerSample() << "b\n";
    std::cout << "\nSamples: "  << wav.getSamplesCount() << "\n";

    frames=wav.getSamplesCount()/SIZE;
    std::cout << "Frames:" << frames << "\n";

    for(k=0; k < frames; k++)
      {
	// Set up frame with zero padding
	for(j=0;j<SIZE;j++)
	  {
	    if(j < SIZE/2)
	      {
		frameLeft[j] = wav.toArray()[j + k*SIZE];
		frameRight[j] = wavR.toArray()[j + k*SIZE];
	      }
	    else 
	      {
		frameLeft[j] = 0;
		frameRight[j] = 0;
	      }
	  }
	
	spectrumLeft = fft->fft(frameLeft);
	spectrumRight = fft->fft(frameRight);    
	
	
	wavRightEnergy = wavLeftEnergy = 0.0;

	for(i=0; i<SIZE; i++) {
	  wavRightEnergy += (spectrumRight[i] * 
			     Aquila::ComplexType(spectrumRight[i].real(), -spectrumRight[i].imag())).real();
	  wavLeftEnergy += (spectrumLeft[i] * 
			    Aquila::ComplexType(spectrumLeft[i].real(), -spectrumLeft[i].imag())).real();
	  correlation[i] = spectrumRight[i]* Aquila::ComplexType(spectrumLeft[i].real(), -spectrumLeft[i].imag());
	}
	
	fft->ifft(correlation, autoCor);
	
#ifdef DEBUG
	std::cout << "\nSamples Left: " <<  spectrumLeft.size();;
	std::cout << "\nSamples Right: " << spectrumRight.size();
	std::cout << "\nEnergy Left: " <<  wavLeftEnergy;
	std::cout << "\nEnergy Right: " << wavRightEnergy;
	std::cout << "\n";
	
	for(int i=0; i<SIZE; i++) {
	  std::cout << "\nwav[" << i << "]=" << frameLeft[i];
	  std::cout << "\nwavR[" << i << "]=" << frameRight[i];
	}
	std::cout << "\n";

	for(int i=0; i<SIZE; i++) {
	  std::cout << "\nCor[" << i << "]=" << autoCor[i];
	}
	std::cout << "\n";
#endif
	
	// Find peak in correlation. Look only within a 1ms window
	maxCorVal=0;
	maxCorLoc=0;
	for(i=0; i <SIZE; i++) {
	  if (autoCor[i] > maxCorVal)
	    {
	      maxCorVal = autoCor[i];
	      maxCorLoc = i;
	    }
	}

	
	// Calculate the time offset
	if (maxCorLoc < SIZE/2.0) { 
	  deltaT = (maxCorLoc/wav.getSampleFrequency()) * 1000.0;
	}
	else 
	  {
	  deltaT = ((maxCorLoc - SIZE)/wav.getSampleFrequency()) * 1000.0;
	  }
	deltaTotal+=deltaT;
        deltaE +=  wavRightEnergy/wavLeftEnergy;


#if 1
	std::cout << "\nDeltaT=" << deltaT << "ms\n" << "DeltaE:" 
		  << (wavLeftEnergy - wavRightEnergy) <<  "SIZE:" << SIZE << "\n";
	

	std::cout << "\nMaxCorLoc= " << maxCorLoc << "\nMax Correlation=" << maxCorVal;
	std::cout << "\nEnergy Left: " <<  wavLeftEnergy;
	std::cout << "\nEnergy Right: " << wavRightEnergy;
 	std::cout << "\ndeltaSamples=" << deltaT*wav.getSampleFrequency()/1000.0 << "\n";
#endif	

	std::cout << "Max Cor Location=" << maxCorLoc << " Value:" << maxCorVal << "\n";
	std::cout << "DeltaT:" << deltaT << "  DeltaT Average=" << deltaTotal/(k+1) << "ms\n";	
	std::cout << "DeltaE Average=" << deltaE/(k+1) << "\n";	
      }
    
    return 0;
}
