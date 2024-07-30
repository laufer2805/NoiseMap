#include <iostream>
#include <iostream>
#include <string.h>
using namespace std;

class Color
{
    private:
        // HexToDec helper function - converts hex string of length 6 to 3 dec numbers saved as color in rgb
        void HexToDec (std::string colorHex)
        {
            //if (hex.length() != 6) return NULL;
            string temp;
            temp = colorHex.substr(0, 2);
            r = stoi(temp, 0, 16);
            temp = colorHex.substr(2, 2);
            g = stoi(temp, 0, 16);
            temp = colorHex.substr(4, 2);
            b = stoi(temp, 0, 16);
        }
    
    public:
        // Values for rgb color need to be ints of range [0, 255]
        int r, g, b;
        int paletteLength;

        // palette from lospec SLSO8, by Luis Miguel Maldonado
        static Color* SLSO8 (void)
        {
            Color* slso8 = new Color[8];
            slso8[0].HexToDec ("0d2b45");
            slso8[1].HexToDec ("203c56");
            slso8[2].HexToDec ("544e68");
            slso8[3].HexToDec ("8d697a");
            slso8[4].HexToDec ("d08159");
            slso8[5].HexToDec ("ffaa5e");
            slso8[6].HexToDec ("ffd4a3");
            slso8[7].HexToDec ("ffecd6");
            slso8[0].paletteLength = 8;
            return slso8;
        }
        
        // palette from lospec Twilight 5 by Star
        static Color* Twilight5 (void)
        {
            Color* twilight5 = new Color[5];
            twilight5[0].HexToDec ("fbbbad");
            twilight5[1].HexToDec ("ee8695");
            twilight5[2].HexToDec ("4a7a96");
            twilight5[3].HexToDec ("333f58");
            twilight5[4].HexToDec ("292831");
            twilight5[0].paletteLength = 5;
            return twilight5;
        }
        
        // palette from lospec OIL 6 by GrafxKid
        static Color* OIL6 (void)
        {
            Color* oil6 = new Color[6];
            oil6[0].HexToDec ("fbf5ef");
            oil6[1].HexToDec ("f2d3ab");
            oil6[2].HexToDec ("c69fa5");
            oil6[3].HexToDec ("8b6d9c");
            oil6[4].HexToDec ("494d7e");
            oil6[5].HexToDec ("272744");
            oil6[0].paletteLength = 6;
            return oil6;
        }
        
        // palette from lospec MOONLIGHT GB by Tofu
        static Color* MoonlightGb (void)
        {
            Color* moonlightGb = new Color[4];
            moonlightGb[0].HexToDec ("0f052d");
            moonlightGb[1].HexToDec ("203671");
            moonlightGb[2].HexToDec ("36868f");
            moonlightGb[3].HexToDec ("5fc75d");
            moonlightGb[0].paletteLength = 4;
            return moonlightGb;
        }

        // palette from lospec chasm by dysphoriaa
        static Color* Chasm (void)
        {
            Color* chasm = new Color[22];
            chasm[0].HexToDec ("85daeb");
            chasm[1].HexToDec ("5fc9e7");
            chasm[2].HexToDec ("5fa1e7");
            chasm[3].HexToDec ("5f6ee7");
            chasm[4].HexToDec ("4c60aa");
            chasm[5].HexToDec ("444774");
            chasm[6].HexToDec ("32313b");
            chasm[7].HexToDec ("463c5e");
            chasm[8].HexToDec ("5d4776");
            chasm[9].HexToDec ("855395");
            chasm[10].HexToDec ("ab58a8");
            chasm[11].HexToDec ("ca60ae");
            chasm[12].HexToDec ("f3a787");
            chasm[13].HexToDec ("f5daa7");
            chasm[14].HexToDec ("8dd894");
            chasm[15].HexToDec ("5dc190");
            chasm[16].HexToDec ("4ab9a3");
            chasm[17].HexToDec ("4593a5");
            chasm[18].HexToDec ("5efdf7");
            chasm[19].HexToDec ("ff5dcc");
            chasm[20].HexToDec ("fdfe89");
            chasm[21].HexToDec ("ffffff");
            chasm[0].paletteLength = 22;
            return chasm;
        }
};

