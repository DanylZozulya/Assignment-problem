#include <algorithm>
#include <iostream>
#include <limits>
#include <list>
#include <string>
#include <vector>
#include <sstream>

using namespace std;


void step1(vector<vector<int>>& matrix, int& step)
    {
        // process rows
        for (auto& row: matrix) {
            auto smallest = *min_element(begin(row), end(row));
            if (smallest > 0)
                for (auto& n: row)
                    n -= smallest;
        }

        // process cols
        int sz = matrix.size(); // square matrix is granted
        for (int j=0; j<sz; ++j) {
            int minval = numeric_limits<int>::max();
            for (int i=0; i<sz; ++i) {
                minval = min(minval, matrix[i][j]);
            }

            if (minval > 0) {
                for (int i=0; i<sz; ++i) {
                    matrix[i][j] -= minval;
                }
            }
        }

        step = 2;
    }

/* helper to clear the temporary vectors */
    inline void clear_covers(vector<int>& cover)
    {
        for (auto& n: cover) n = 0;
    }

    void step2(const vector<vector<int>>& matrix,
               vector<vector<int>>& M,
               vector<int>& RowCover,
               vector<int>& ColCover,
               int& step)
    {
        int sz = matrix.size();

        for (int r=0; r<sz; ++r)
            for (int c=0; c<sz; ++c)
                if (matrix[r][c] == 0)
                    if (RowCover[r] == 0 && ColCover[c] == 0) {
                        M[r][c] = 1;
                        RowCover[r] = 1;
                        ColCover[c] = 1;
                    }

        clear_covers(RowCover); // reset vectors for posterior using
        clear_covers(ColCover);

        step = 3;
    }


    void step3(const vector<vector<int>>& M,
               vector<int>& ColCover,
               int& step)
    {
        int sz = M.size();
        int colcount = 0;

        for (int r=0; r<sz; ++r)
            for (int c=0; c<sz; ++c)
                if (M[r][c] == 1)
                    ColCover[c] = 1;

        for (auto& n: ColCover)
            if (n == 1)
                colcount++;

        if (colcount >= sz) {
            step = 7; // solution found
        }
        else {
            step = 4;
        }
    }

// Following functions to support step 4
    void find_a_zero(int& row,
                     int& col,
                     const vector<vector<int>>& matrix,
                     const vector<int>& RowCover,
                     const vector<int>& ColCover)
    {
        int r = 0;
        int c = 0;
        int sz = matrix.size();
        bool done = false;
        row = -1;
        col = -1;

        while (!done) {
            c = 0;
            while (true) {
                if (matrix[r][c] == 0 && RowCover[r] == 0 && ColCover[c] == 0) {
                    row = r;
                    col = c;
                    done = true;
                }
                c += 1;
                if (c >= sz || done)
                    break;
            }
            r += 1;
            if (r >= sz)
                done = true;
        }
    }

    bool star_in_row(int row,
                     const vector<vector<int>>& M)
    {
        bool tmp = false;
        for (unsigned c = 0; c < M.size(); c++)
            if (M[row][c] == 1)
                tmp = true;

        return tmp;
    }


    void find_star_in_row(int row,
                          int& col,
                          const vector<vector<int>>& M)
    {
        col = -1;
        for (unsigned c = 0; c < M.size(); c++)
            if (M[row][c] == 1)
                col = c;
    }


    void step4(const vector<vector<int>>& matrix,
               vector<vector<int>>& M,
               vector<int>& RowCover,
               vector<int>& ColCover,
               int& path_row_0,
               int& path_col_0,
               int& step)
    {
        int row = -1;
        int col = -1;
        bool done = false;

        while (!done){
            find_a_zero(row, col, matrix, RowCover, ColCover);

            if (row == -1){
                done = true;
                step = 6;
            }
            else {
                M[row][col] = 2;
                if (star_in_row(row, M)) {
                    find_star_in_row(row, col, M);
                    RowCover[row] = 1;
                    ColCover[col] = 0;
                }
                else {
                    done = true;
                    step = 5;
                    path_row_0 = row;
                    path_col_0 = col;
                }
            }
        }
    }

// Following functions to support step 5
    void find_star_in_col(int c,
                          int& r,
                          const vector<vector<int>>& M)
    {
        r = -1;
        for (unsigned i = 0; i < M.size(); i++)
            if (M[i][c] == 1)
                r = i;
    }

    void find_prime_in_row(int r,
                           int& c,
                           const vector<vector<int>>& M)
    {
        for (unsigned j = 0; j < M.size(); j++)
            if (M[r][j] == 2)
                c = j;
    }

    void augment_path(vector<vector<int>>& path,
                      int path_count,
                      vector<vector<int>>& M)
    {
        for (int p = 0; p < path_count; p++)
            if (M[path[p][0]][path[p][1]] == 1)
                M[path[p][0]][path[p][1]] = 0;
            else
                M[path[p][0]][path[p][1]] = 1;
    }

    void erase_primes(vector<vector<int>>& M)
    {
        for (auto& row: M)
            for (auto& val: row)
                if (val == 2)
                    val = 0;
    }


    void step5(vector<vector<int>>& path,
               int path_row_0,
               int path_col_0,
               vector<vector<int>>& M,
               vector<int>& RowCover,
               vector<int>& ColCover,
               int& step)
    {
        int r = -1;
        int c = -1;
        int path_count = 1;

        path[path_count - 1][0] = path_row_0;
        path[path_count - 1][1] = path_col_0;

        bool done = false;
        while (!done) {
            find_star_in_col(path[path_count - 1][1], r, M);
            if (r > -1) {
                path_count += 1;
                path[path_count - 1][0] = r;
                path[path_count - 1][1] = path[path_count - 2][1];
            }
            else {done = true;}

            if (!done) {
                find_prime_in_row(path[path_count - 1][0], c, M);
                path_count += 1;
                path[path_count - 1][0] = path[path_count - 2][0];
                path[path_count - 1][1] = c;
            }
        }

        augment_path(path, path_count, M);
        clear_covers(RowCover);
        clear_covers(ColCover);
        erase_primes(M);

        step = 3;
    }

// methods to support step 6
    void find_smallest(int& minval,
                       const vector<vector<int>>& matrix,
                       const vector<int>& RowCover,
                       const vector<int>& ColCover)
    {
        for (unsigned r = 0; r < matrix.size(); r++)
            for (unsigned c = 0; c < matrix.size(); c++)
                if (RowCover[r] == 0 && ColCover[c] == 0)
                    if (minval > matrix[r][c])
                        minval = matrix[r][c];
    }

    void step6(vector<vector<int>>& matrix,
               const vector<int>& RowCover,
               const vector<int>& ColCover,
               int& step)
    {
        int minval = numeric_limits<int>::max();
        find_smallest(minval, matrix, RowCover, ColCover);

        int sz = matrix.size();
        for (int r = 0; r < sz; r++)
            for (int c = 0; c < sz; c++) {
                if (RowCover[r] == 1)
                    matrix[r][c] += minval;
                if (ColCover[c] == 0)
                    matrix[r][c] -= minval;
            }

        step = 4;
    }


vector<vector<int>> takeMatrix(string strArr[], int arrLength){
    vector<vector<int>> matrix(arrLength);
    string num;
    int temp;
    for (int x = 0; x < arrLength; x++)
    {
        for (int y = 0; y < strArr[x].length(); y++)
        {
            if ((strArr[x][y] == '(' || strArr[x][y] == ')' || strArr[x][y] == ',') && !num.empty())
            {
                istringstream(num) >> temp;
                matrix[x].push_back(temp);
                num.clear();
            }
            else if (isdigit(strArr[x][y]))
            {
                num += strArr[x][y];
            }
        }
    }
    return matrix;
}

string OptimalAssignments(string strArr[], int arrLength) {
    vector<vector<int>> matrix = takeMatrix(strArr, arrLength);
    vector<vector<int>> M (arrLength, vector<int>(arrLength, 0));

    /*Define two vectors RowCover and ColCover that are used to "cover"
     *the rows and columns of the cost matrix */
    vector<int> RowCover (arrLength, 0);
    vector<int> ColCover (arrLength, 0);

    int path_row_0, path_col_0; //temporary to hold the smallest uncovered value

    // Array for the augmenting path algorithm
    vector<vector<int>> path (arrLength+1, vector<int>(2, 0));

    bool done = false;
    int step = 1;
    while (!done) {
        switch (step) {
            case 1:
                step1(matrix, step);
                break;
            case 2:
                step2(matrix, M, RowCover, ColCover, step);
                break;
            case 3:
                step3(M, ColCover, step);
                break;
            case 4:
                step4(matrix, M, RowCover, ColCover, path_row_0, path_col_0, step);
                break;
            case 5:
                step5(path, path_row_0, path_col_0, M, RowCover, ColCover, step);
                break;
            case 6:
                step6(matrix, RowCover, ColCover, step);
                break;
            case 7:
                for (auto& vec: M) {vec.resize(matrix.begin()->size());}
                M.resize(arrLength);
                done = true;
                break;
            default:
                done = true;
                break;
        }
    }

    string res;

    for (int i = 0;i < arrLength;i++){
        for (int j = 0;j < arrLength;j++){
            if (M[i][j] == 1){
                res += "(" + to_string(i+1) + "-" + to_string(j+1) + ")";
            }
        }
    }
    return res;


}

template<typename C>
auto sum(const C& cnt)
{
    auto it = cnt.begin();
    auto ret = *it;
    for(++it;it!=cnt.end();it++)
        ret += static_cast<int>(*it);
    return static_cast<int>(ret);
}

int main() {
    sstr
    return 0;
}
