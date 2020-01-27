#include <iostream>
#include <math.h>

using namespace std;

// Read this paper : https://www.st.com/content/ccc/resource/technical/document/application_note/a0/f0/a0/62/3b/69/47/66/DM00119044.pdf/files/DM00119044.pdf/jcr:content/translations/en.DM00119044.pdf

float Y[6][3] = {
    {0,0,1},    //ideal Z axis up
    {0,0,-1},   //ideal Z axis down
    {0,1,0},    //ideal Y axis up
    {0,-1,0},   //ideal Y axis down
    {1,0,0},    //ideal X axis up
    {-1,0,0},   //ideal X axis down
};

float w[6][4] = {
    {0.03, -0.02, 0.96, 1.0},   // measured Z axis up
    {0.03, 0.02, -1.04, 1.0},   // measured Z axis down
    {0.03, 0.97, 0.01, 1.0},    // measured Y axis up
    {-0.03, -1.03, -0.01, 1.0}, // measured Y axis down
    {0.99, -0.03, -0., 1.0},    // measured X axis up
    {-1.01, -0.03, 0., 1.0}     // measured X axis down
};

bool invertMatrix(const float m[16], float invOut[16])
{
    float inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;
    return true;
}

int main()
{
    float w_T[4][6];
    for(int row=0; row < 6; row++){
        for(int col=0; col < 4; col++){
            w_T[col][row]=0;
            w_T[col][row]=w[row][col];
        }
    }

    float dot_prod[4][4];
    for(int r = 0; r < 4; ++r)
        for(int c = 0; c < 4; ++c){
            dot_prod[r][c] = 0.0;
            for(int k = 0; k < 6; ++k){
                dot_prod[r][c] += w_T[r][k] * w[k][c];
            }
        }

    float reshaped[16];
    float inverted[16];
    int idx=0;
    for(int r=0; r < 4; r++){
        for(int c=0; c < 4; c++){
            reshaped[idx] = 0;
            inverted[idx] = 0;
            reshaped[idx++] = dot_prod[r][c];
        }
    }

    invertMatrix(reshaped, inverted);

    idx=0;
    float inverted_reshaped[4][4];
    for(int r=0; r < 4; r++){
        for(int c=0; c < 4; c++){
            inverted_reshaped[r][c] = 0.;
            inverted_reshaped[r][c] = inverted[idx++];
        }
    }

    float dot_prod2[4][6];
    for(int r = 0; r < 4; ++r)
        for(int c = 0; c < 6; ++c){
            dot_prod2[r][c] = 0.0;
            for(int k = 0; k < 4; ++k){
                dot_prod2[r][c] += inverted_reshaped[r][k] * w_T[k][c];
            }
        }

    float X[4][3];
    for(int r = 0; r < 4; ++r)
        for(int c = 0; c < 3; ++c){
            X[r][c] = 0.0;
            for(int k = 0; k < 6; ++k){
                X[r][c] += dot_prod2[r][k] * Y[k][c];
            }
        }


    return 0;
}
