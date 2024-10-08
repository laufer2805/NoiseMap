#include <iostream>
#include <fstream>
#include <string>
#include </media/data/programming/GitHub/NoiseMap/Scripts/PseudoRandomNumberGenerator.cpp>
#include </media/data/programming/GitHub/NoiseMap/Scripts/Color.cpp>

class NoiseMap
{
    private:
        Color* palette;
        int paletteLength;
        int* seedWordInt;
        int vowelCount, consonantCount, seedWordLength, wordCount;
        double averageLetter;
        PseudoRandomNumberGenerator prng;

        void ConvertSeedWordToInt ()
        {
            seedWordInt = new int[seedWordLength];
            for (int i = 0; i < seedWordLength; i++)
            {
                seedWordInt[i] = (int)seedWord[i];
            }
            SeedWordAnalysis ();
        }

        void SeedWordAnalysis()
        {
            seedWordLength = seedWord.length();
            vowelCount = 0;
            consonantCount = 0;
            wordCount = 1;
            averageLetter = 0; 
            for (int i = 0; i < seedWordLength; i++)
            {
                if (seedWord[i] == 'A' || seedWord[i] == 'a' || seedWord[i] == 'E' || seedWord[i] == 'e' || seedWord[i] == 'I' || seedWord[i] == 'i' || seedWord[i] == 'O' || seedWord[i] == 'o' || seedWord[i] == 'U' || seedWord[i] == 'u' )
                {
                    vowelCount++;
                }
                if (seedWord[i] == ' ')
                {
                    wordCount++;
                }
                averageLetter += seedWord[i];
            }
            averageLetter /= (double)seedWordLength;
            consonantCount = seedWordLength - vowelCount;
        }
        
        // draw functions: DrawCircle (int centerX, int centerY, int radius, bool isFilled)
        /*
            draw functions: 
                DrawCircle (int centerX, int centerY, int radius, bool isFilled)
                DrawLine (int startX, int startY, int endX, int endY)
                
        */
        void DrawLine (int startX, int startY, int endX, int endY, int colorIndex)
        {
            for (int i = startX; i <= endX; i++)
            {
                int j = (endY - startY)/(endX - startX) * (i - endX) + endY;
                image[i][j] = colorIndex;
            }
        }
        void DrawCircle (int centerX, int centerY, int radius, bool isFilled, int colorIndex)
        {
            for (int i = 0; i < sizeX; i++)
            {
                for (int j = 0; j < sizeY; j++)
                {
                    int deltaX = i - centerX, deltaY = j - centerY;
                    if (deltaX * deltaX + deltaY * deltaY == radius * radius)
                    {
                        image[i][j] = colorIndex;
                    }
                    if (isFilled == true)
                    {
                        if (deltaX * deltaX + deltaY * deltaY < radius * radius)
                        {
                            image[i][j] = colorIndex;
                        }
                    }
                }
            }
        }

        void DrawRandomLines (int numberOfLines, int colorIndex)
        {
            double *startsX, *startsY, *endsX, *endsY;
            // need to check if starts and ends are good
            double *tempSeq;
            tempSeq = prng.GenerateSequence2D (averageLetter, (double)wordCount, 8);
            startsX = prng.GenerateSequence2D (tempSeq[0], tempSeq[1], numberOfLines);
            startsY = prng.GenerateSequence2D (tempSeq[2], tempSeq[3], numberOfLines);
            endsX = prng.GenerateSequence2D (tempSeq[4], tempSeq[5], numberOfLines);
            endsY = prng.GenerateSequence2D (tempSeq[6], tempSeq[7], numberOfLines);
            for (int i = 0; i < numberOfLines; i++)
            {
                startsX[i] = Clamp (startsX[i], 0, sizeX);
                endsX[i] = Clamp (endsX[i], 0, sizeX);
                startsY[i] = Clamp (startsY[i], 0, sizeY);
                endsY[i] = Clamp (endsY[i], 0, sizeY);
                DrawLine (startsX[i], startsY[i], endsX[i], endsY[i], colorIndex);
            }
        }

        void DrawRandomCircles (int numberOfCircles, int colorIndex)
        {
            // needed random numbers: centerX, centerY, radius
            /*
                algo:
                radii depend on numberOfCircles - min = 5, max = minSize / numberOfCircles
            */
            double *centersX, *centersY, *radii;
            double *tempSeq;
            tempSeq = prng.GenerateSequence2D ((double)vowelCount, (double)wordCount, 6);
            centersX = prng.GenerateSequence2D(tempSeq[0], tempSeq[1], numberOfCircles);
            centersY = prng.GenerateSequence2D(tempSeq[2], tempSeq[3], numberOfCircles);
            radii = prng.GenerateSequence2D(tempSeq[4], tempSeq[5], numberOfCircles);
            double minSize = Min((double)sizeX, (double)sizeY);
            for (int i = 0; i < numberOfCircles; i++)
            {
                centersX[i] = Clamp (centersX[i], 0, sizeX);
                centersY[i] = Clamp (centersY[i], 0, sizeY);
                radii[i] = Clamp (radii[i], 5, minSize / (double)numberOfCircles);
                DrawCircle ((int)centersX[i], (int)centersY[i], (int)radii[i], true, colorIndex);
            }
        }

        void GenerateArray ()
        {
            image = new int*[sizeX];
            for (int i = 0; i < sizeX; i++)
            {
                image[i] = new int[sizeY];
            }
        }
        bool PaletteChooser (int paletteName)
        {
            if (paletteName == 1)
            {
                palette = Color::SLSO8();
            }
            else if (paletteName == 2)
            {
                palette = Color::Twilight5 ();
            }
            else if (paletteName == 3)
            {
                palette = Color::OIL6 ();
            }
            else if (paletteName == 4)
            {
                palette = Color::MoonlightGb ();
            }
            else if (paletteName == 5)
            {
                palette = Color::Chasm ();
            }
            else
            {
                std::cout << "Wrong palette chosen!" << std::endl;
                return false;
            }
            paletteLength = palette[0].paletteLength;
            return true;
        }
        void InitialFill ()
        {
            for (int i = 0; i < sizeX; i++)
            {
                for (int j = 0; j < sizeY; j++)
                {
                    image[i][j] = prng.IntRange2D (i, j, 0, paletteLength);
                }
            }
        }
        void PrintArray (void)
        {
            for (int i = 0; i < sizeX; i++)
            {
                for (int j = 0; j < sizeY; j++)
                {
                    std::cout << image[i][j] << "\t";
                }
                std::cout << std::endl;
            }
        }
        void PrintToFile (string name)
        {
            std::ofstream outFile;
            // need to check this part of code - should work on c++ 17
            /*
                std::filesystem::path cwd = std::filesystem::current_path();
                outFile.open(cwd + name + ".ppm");
            */
            outFile.open("/media/data/programming/GitHub/NoiseMap/Archive/" + name + ".ppm");
            outFile << "P3" << endl << sizeX << " " << sizeY << endl << 255 << endl;
            for (int i = 0; i < sizeX; i++)
            {
                for (int j = 0;j < sizeY; j++)
                {
                    outFile << palette[image[i][j]].r << " " << palette[image[i][j]].g << " " << palette[image[i][j]].b << endl;
                }
            }
        }
    public:
        int sizeX, sizeY;           // size of image in pixels
        int** image;
        std::string seedWord;

        void GenerateNoiseMap (void)
        {
            // inputs
            std::cout << "Input seed word: ";
            std::getline(std::cin, seedWord);
            std::cout << "Input sizeX: ";
            std::cin >> sizeX;
            std::cout << "Input sizeY: ";
            std::cin >> sizeY;
            std::cout << "Choose palette:\n(1) SLSO8\n(2) Twilight5\n(3) OIL6\n(4) MoonlightGb\n(5) Chasm\n";
            bool flag = false;
            int paletteName;
            while (true)
            {
                cin >> paletteName;
                PaletteChooser (paletteName);
                std::cout << "paletteLength = " << paletteLength;
                if (paletteName > 0 && paletteName < 6)
                {
                    break;
                }
            }

            // generating image
            GenerateArray ();
            std::cout << "Array Generated" << std::endl;
            ConvertSeedWordToInt ();
            InitialFill ();
            std::cout << "Initial fill done" << std::endl;
            DrawRandomCircles (vowelCount, consonantCount % paletteLength);
            DrawRandomLines ((int)averageLetter, wordCount % paletteLength);
            std::cout << "DrawRandomCircles done" << std::endl;
            
            // output
            PrintToFile (seedWord);
            std::cout << "File generated" << std::endl;
        }
        
};

int main (void)
{
    NoiseMap nm;
    nm.GenerateNoiseMap ();
}