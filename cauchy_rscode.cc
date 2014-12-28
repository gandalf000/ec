/**
 * Copyright (c) 2014, The Authors. All rights reserved.
 * @file cauchy_rscode.cpp
 * @brief Cauchy Reed-Solomon encoding and decoding library
 * @author zhouw(zhouw@zhigu.com)
 * @date 2014-12-29
 */

#include "common/cauchy_rscode.h"
#include <stdint.h>
#include <string.h>

namespace phenix {
namespace common {

static inline void MemXor(const char* src1,  // source buffer 1
                        const char *src2,  // source buffer 2
                        char *dst,         // Xor src1 and src2 (dst = src1 ^ src2)
                        int size) {        // buffer size
    const int64_t *l1 = (const int64_t *)src1;
    const int64_t *l2 = (const int64_t *)src2;
    int64_t *l3 = reinterpret_cast<int64_t *>(dst);

    for (int count = 0; count < size; count += sizeof(l3[0]), l1++, l2++, l3++) {
        *l3 = ((*l1)  ^ (*l2));
    }
}

void CauchyRSCoder::_FreeSchedule(int **schedule) {
    int i = 0;
    for ( ; schedule[i][0] >= 0; i++)
        delete[] schedule[i];
    delete[] schedule[i];
    delete[] schedule;
}

int CauchyRSCoder::_CountCauchyOnes(int num) {
    int ones_count = 0;
    static int PPs = -1;
    static int ONEs[kWordBits];
    static int NOs = 0;

    int highbit = (1 << (kWordBits - 1));
    if (PPs == -1) {
        ones_count = 0;
        PPs = m_galois_operator->Multiply(highbit, 2);
        for (int i = 0; i < kWordBits; i++) {
            if (PPs & (1 << i)) {
                ONEs[ones_count] = (1 << i);
                ones_count++;
            }
        }
        NOs = ones_count;
    }

    ones_count = 0;
    for (int i = 0; i < kWordBits; i++)
        if (num & (1 << i))
            ones_count++;

    int cur_ones_count = 0;
    cur_ones_count = ones_count;
    for (int i = 1; i < kWordBits; i++) {
        if (num & highbit) {
            num ^= highbit;
            num <<= 1;
            num ^= PPs;
            cur_ones_count--;
            for (int j = 0; j < NOs; j++) {
                cur_ones_count += (num & ONEs[j]) ? 1 : -1;
            }
        } else {
          num <<= 1;
        }
        ones_count += cur_ones_count;
    }
    return ones_count;
}

int *CauchyRSCoder::_GenerateEncodeMatrix() {
    // generate original cauchy coding matrix
    int *matrix = new int[sizeof(int) * m_num_data_parts * m_num_code_parts]; // NOLINT
    int index = 0, tmp = 0;
    for (int i = 0; i < m_num_code_parts; i++) {
        for (int j = 0; j < m_num_data_parts; j++) {
            index = i * m_num_data_parts + j;
            matrix[index] = m_galois_operator->Divide(1, (i ^ (m_num_code_parts + j)));
        }
    }

    // improve the cauchy coding matrix, make the first line of coding matrix all one
    for (int i = 0; i < m_num_data_parts; i++) {
        if (matrix[i] != 1) {
            tmp = m_galois_operator->Divide(1, matrix[i]);
            index = i;
            for (int j = 0; j < m_num_code_parts; j++) {
                matrix[index] = m_galois_operator->Multiply(matrix[index], tmp);
                index += m_num_data_parts;
            }
        }
    }

    // improve other lines of coding matrix, we do the following:
    // 1 Count the number of ones in the bit representation of the row
    // 2 Count the number of ones in the bit representation of the row
    //   divided by element M[i, j] for each j
    // 3 Whichever value of j gives the minimal number of ones, if it
    //   improves the number of ones in the original row, divide row i by M[i, j].
    int min_ones_count = 0;
    int cur_ones_count = 0;
    int best_m_index = -1;
    for (int i = 1; i < m_num_code_parts; i++) {
        // tep 1
        min_ones_count = 0;
        index = i * m_num_data_parts;
        for (int j = 0; j < m_num_data_parts; j++)
            min_ones_count += _CountCauchyOnes(matrix[index + j]);

        // tep 2
        best_m_index = -1;
        for (int j = 0; j < m_num_data_parts; j++) {
            if (matrix[index + j] != 1) {
                int tmp = m_galois_operator->Divide(1, matrix[index + j]);
                cur_ones_count = 0;
                for (int k = 0; k < m_num_data_parts; k++) {
                    cur_ones_count += _CountCauchyOnes(m_galois_operator->Multiply(
                                                        matrix[index + k], tmp));
                }

                if (cur_ones_count < min_ones_count) {
                    min_ones_count = cur_ones_count;
                    best_m_index = j;
                }
            }
        }

        // tep 3
        if (best_m_index != -1) {
            tmp = m_galois_operator->Divide(1, matrix[index + best_m_index]);
            for (int j = 0; j < m_num_data_parts; j++)
                matrix[index + j] = m_galois_operator->Multiply(matrix[index + j], tmp);
        }
    }
    return matrix;
}

char *CauchyRSCoder::_MatrixToBitMatrix(int *matrix) {
    size_t matrix_size = m_num_data_parts * m_num_code_parts * kWordBits * kWordBits;
    char *bit_matrix = new char[matrix_size];
    int num_bits_per_row = m_num_data_parts * kWordBits;
    int matrix_element = 0;
    int bit_matrix_row = 0;
    int bit_matrix_col = 0;
    for (int i = 0; i < m_num_code_parts; i++) {
        bit_matrix_col = bit_matrix_row;
        for (int j = 0; j < m_num_data_parts; j++) {
            matrix_element = matrix[i * m_num_data_parts + j];
            // xpand every element of matrix to w * w bit matrix
            for (int m = 0; m < kWordBits; m++) {
                for (int n = 0; n < kWordBits; n++) {
                    bit_matrix[bit_matrix_col + m + n * num_bits_per_row]
                        = ((matrix_element & (1 << n)) ? 1 : 0);
                }

                matrix_element = m_galois_operator->Multiply(matrix_element, 2);
            }
            bit_matrix_col += kWordBits;
        }
        bit_matrix_row += num_bits_per_row * kWordBits;
    }

    return bit_matrix;
}

int **CauchyRSCoder::_BitMatrixToSchedule(int num_data_parts,
                                        int num_code_parts,
                                        char *bit_matrix) {
    size_t size = num_code_parts * kWordBits;
    int *diff = new int[size];
    int *from = new int[size];
    int *flink = new int[size];
    int *blink = new int[size];

    int num_ones = 0, num_ops = 0;
    int best_diff = num_data_parts * kWordBits + 1;
    int best_row_index = -1;
    char *matrix_iter = bit_matrix;
    for (int i = 0; i < num_code_parts * kWordBits; i++) {
        num_ones = 0;
        for (int j = 0; j < num_data_parts * kWordBits; j++) {
            if (*matrix_iter++ == 1) {
                num_ones++;
            }
        }

        diff[i] = num_ones;
        from[i] = -1;
        flink[i] = i + 1;
        blink[i] = i - 1;
        if (num_ones < best_diff) {
            best_diff = num_ones;
            best_row_index = i;
        }
    }
    flink[num_code_parts * kWordBits - 1] = -1;

    size = num_data_parts * num_code_parts * kWordBits * kWordBits + 1;
    int **operations = new int*[size];
    int top = 0, row_index = 0, op_to_do = 0;
    while (top != -1) {
        row_index = best_row_index;

        if (blink[row_index] == -1) {
            top = flink[row_index];
            if (top != -1) {
                blink[top] = -1;
            }
        } else {
            flink[blink[row_index]] = flink[row_index];
            if (flink[row_index] != -1) {
                blink[flink[row_index]] = blink[row_index];
            }
        }

        matrix_iter = bit_matrix + row_index * num_data_parts * kWordBits;
        if (from[row_index] == -1) {
            op_to_do = 0;
            for (int j = 0; j < num_data_parts * kWordBits; j++) {
                if (matrix_iter[j]) {
                    operations[num_ops] = new int[5];
                    operations[num_ops][4] = op_to_do;
                    operations[num_ops][0] = j / kWordBits;
                    operations[num_ops][1] = j % kWordBits;
                    operations[num_ops][2] = num_data_parts + row_index / kWordBits;
                    operations[num_ops][3] = row_index % kWordBits;
                    op_to_do = 1;
                    num_ops++;
                }
            }
        } else {
            operations[num_ops] = new int[5];
            operations[num_ops][4] = 0;
            operations[num_ops][0] = num_data_parts + from[row_index]/kWordBits;
            operations[num_ops][1] = from[row_index] % kWordBits;
            operations[num_ops][2] = num_data_parts + row_index/kWordBits;
            operations[num_ops][3] = row_index % kWordBits;
            num_ops++;
            char *b1 = bit_matrix + from[row_index] * num_data_parts * kWordBits;
            for (int j = 0; j < num_data_parts*kWordBits; j++) {
                if (matrix_iter[j] ^ b1[j]) {
                    operations[num_ops] = new int[5];
                    operations[num_ops][4] = 1;
                    operations[num_ops][0] = j / kWordBits;
                    operations[num_ops][1] = j % kWordBits;
                    operations[num_ops][2] = num_data_parts + row_index/kWordBits;
                    operations[num_ops][3] = row_index % kWordBits;
                    op_to_do = 1;
                    num_ops++;
                }
            }
        }
        best_diff = num_data_parts*kWordBits+1;
        for (int i = top; i != -1; i = flink[i]) {
            num_ones = 1;
            char *b1 = bit_matrix + i*num_data_parts*kWordBits;
            for (int j = 0; j < num_data_parts * kWordBits; j++)
                num_ones += (matrix_iter[j] ^ b1[j]);
            if (num_ones < diff[i]) {
                from[i] = row_index;
                diff[i] = num_ones;
            }
            if (diff[i] < best_diff) {
                best_diff = diff[i];
                best_row_index = i;
            }
        }
    }

    operations[num_ops] = new int[5];
    operations[num_ops][0] = -1;
    delete[] from;
    delete[] diff;
    delete[] blink;
    delete[] flink;

    return operations;
}

void CauchyRSCoder::_DoScheduleOperations(int **schedule, char **ptrs, int size) {
    int unit_size = kPacketSize * kWordBits;
    for (int count = 0; count < size; count += unit_size) {
        for (int i = 0; schedule[i][0] >= 0; i++) {
            char *src = ptrs[schedule[i][0]] + schedule[i][1] * kPacketSize;
            char *dst = ptrs[schedule[i][2]] + schedule[i][3] * kPacketSize;
            if (schedule[i][4]) {
                MemXor(src, dst, dst, kPacketSize);
            } else {
                memcpy(dst, src, kPacketSize);
            }
        }
        for (int i = 0; i < m_num_data_parts + m_num_code_parts; i++)
            ptrs[i] += unit_size;
    }
}

void CauchyRSCoder::_Init() {
    // generate coding matrix and make it sparse
    int *coding_matrix = _GenerateEncodeMatrix();

    // convert matrix to bitmatrix, to convert multply and divide opertaion on
    // GF(2^8) to more efficient XOR opertion, thus making encoding and decoding
    // much faster
    m_encoding_bit_matrix = _MatrixToBitMatrix(coding_matrix);

    // convert bitmatrix to schedule, to avoid traversing the matrix during encoding
    // schedule is a list of 5-tuples: < op, sd, sb, dd, db >,where op is an operation
    // code: 0 for copy and 1 for XOR, sd is the id of the source device and sb is the bit
    // of the source device. The last two elements, dd and db are the destination device and bit.
    m_encoding_schedule = _BitMatrixToSchedule(m_num_data_parts,
                                                m_num_code_parts, m_encoding_bit_matrix);

    delete[] coding_matrix;
}

void CauchyRSCoder::Encode(char **data_ptrs, char **coding_ptrs, int size) {
    assert(size > 0);
    assert(data_ptrs != NULL);
    assert(coding_ptrs != NULL);

    char *ptrs[m_num_data_parts + m_num_code_parts];
    for (int i = 0; i < m_num_data_parts; i++) {
        ptrs[i] = const_cast<char*>(data_ptrs[i]);
    }

    for (int i = 0; i < m_num_code_parts; i++) {
        ptrs[i + m_num_data_parts] = coding_ptrs[i];
    }
    // do encoding
    _DoScheduleOperations(m_encoding_schedule, ptrs, size);
}

static inline void _InvertBitMatrix(char *matrix, char *inverse, int num_rows) {
    int num_cols = num_rows;

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            inverse[i * num_cols + j] = (i == j) ? 1 : 0;
        }
    }

    /* first -- convert into upper triangular */
    for (int i = 0; i < num_cols; i++) {
        /* swap num_rows if we have a zero i,i element.  if we can't swap, then the
           matrixrix was not inverseertible */

        if ((matrix[i*num_cols + i]) == 0) {
            int j = 0;
            for (j = i + 1; j < num_rows && (matrix[j * num_cols + i]) == 0; j++) {
            }

            for (int k = 0; k < num_cols; k++) {
                char tmp = matrix[i * num_cols + k];
                matrix[i * num_cols + k] = matrix[j * num_cols + k];
                matrix[j * num_cols+k] = tmp;
                tmp = inverse[i * num_cols + k];
                inverse[i * num_cols + k] = inverse[j * num_cols+k];
                inverse[j * num_cols + k] = tmp;
            }
        }

        /* now for each j>i, add a_ji*ai to aj */
        for (int j = i + 1; j != num_rows; j++) {
            if (matrix[j * num_cols + i] != 0) {
                for (int k = 0; k < num_cols; k++) {
                    matrix[j*num_cols+k] ^= matrix[i*num_cols+k];
                    inverse[j*num_cols+k] ^= inverse[i*num_cols+k];
                }
            }
        }
    }

    /* now the matrixrix is upper triangular.  start at the top and multiply down */

    for (int i = num_rows - 1; i >= 0; i--) {
        for (int j = 0; j < i; j++) {
            if (matrix[j * num_cols + i]) {
                for (int k = 0; k < num_cols; k++) {
                    matrix[j * num_cols + k] ^= matrix[i * num_cols+k];
                    inverse[j * num_cols + k] ^= inverse[i * num_cols+k];
                }
            }
        }
    }
}

void CauchyRSCoder::Decode(bool *erased,
                            char **data_ptrs,
                            char **coding_ptrs,
                            int size) {
    assert(size > 0);
    assert(data_ptrs != NULL);
    assert(coding_ptrs != NULL);
    assert(erased != NULL);

    int num_total_parts = m_num_data_parts + m_num_code_parts;
    int good_parts_count = 0;
    for (int i = 0; i < num_total_parts; i++) {
        if (!erased[i]) {
           good_parts_count++;
        }
    }

    assert(good_parts_count >= m_num_data_parts);

    // f there is no erased parts, do nothing
    if (good_parts_count == m_num_data_parts + m_num_code_parts) {
        return;
    }


    /* Preapre, set up ptrs.  It will be as follows:

       - If data drive i has not eraseded, then ptrs[i] = data_ptrs[i].
       - If data drive i has eraseded, then ptrs[i] = coding_ptrs[j], where j is the
            lowest unused non-eraseded coding drive.
       - Elements num_data_parts to num_data_parts+num_erased_data_pars-1 are data_ptrs[]
         of the eraseded data drives.
       - Elements num_data_parts+num_erased_data_pars to num_data_parts+num_erased_data_pars
         +num_erased_code_parts-1 are coding_ptrs[] of the eraseded data drives.

       The array rowid_to_partidx used to map matrix row id to part index;
       The array partidx_to_rowid used to map part index to matrix row id;
     */
    char *ptrs[num_total_parts];
    int rowid_to_partidx[num_total_parts]; // NOLINT
    int partidx_to_rowid[num_total_parts]; // NOLINT
    int good_code_part_index = m_num_data_parts;
    int erased_part_index = m_num_data_parts;
    int num_erased_data_parts = 0;
    int num_erased_code_parts = 0;

    for (int i = 0; i < m_num_data_parts; i++) {
        if (!erased[i]) {
            ptrs[i] = data_ptrs[i];
            rowid_to_partidx[i] = i;
            partidx_to_rowid[i] = i;
        } else {
            while (erased[good_code_part_index]) good_code_part_index++;
            ptrs[i] = coding_ptrs[good_code_part_index - m_num_data_parts];
            rowid_to_partidx[i] = good_code_part_index;
            partidx_to_rowid[good_code_part_index] = i;
            good_code_part_index++;

            ptrs[erased_part_index] = data_ptrs[i];
            rowid_to_partidx[erased_part_index] = i;
            partidx_to_rowid[i] = erased_part_index;
            erased_part_index++;
            num_erased_data_parts++;
        }
    }
    for (int i = m_num_data_parts; i < num_total_parts; i++) {
        if (erased[i]) {
            ptrs[erased_part_index] = coding_ptrs[i - m_num_data_parts];
            rowid_to_partidx[erased_part_index] = i;
            partidx_to_rowid[i] = erased_part_index;
            erased_part_index++;
            num_erased_code_parts++;
        }
    }

    /* Now, we're going to create one decoding matrix which is going to
     * decode erased parts. This matrix has kWordBits * kWordBits
     * * (num_erased_data+num_erased_code) * num_data_parts rows
     */

    char *decoding_bit_matrix = new char[m_num_data_parts * kWordBits * kWordBits
                                    * (num_erased_data_parts + num_erased_code_parts)];

    /* First, if any data drives have eraseded, then initialize the first
     * num_erased_data_parts*kWordBits rows of the decoding matrix from the
     * standard decoding * matrix inversion */
    if (num_erased_data_parts > 0) {
        char *tmp_bit_matrix = new char[m_num_data_parts * m_num_data_parts
                                        * kWordBits * kWordBits];
        char *iter = tmp_bit_matrix;
        for (int i = 0; i < m_num_data_parts; i++) {
            if (rowid_to_partidx[i] == i) {
                // if the corresponding part is not erased, make diagonal matrix
                bzero(iter, m_num_data_parts * kWordBits * kWordBits * sizeof(char)); // NOLINT
                for (int j = 0; j < kWordBits; j++) {
                    iter[j + i*kWordBits + j*m_num_data_parts*kWordBits] = 1;
                }
            } else {
                memcpy(iter, m_encoding_bit_matrix + m_num_data_parts * kWordBits
                        * kWordBits * (rowid_to_partidx[i] - m_num_data_parts),
                        m_num_data_parts * kWordBits * kWordBits * sizeof(char)); // NOLINT
            }
            iter += (m_num_data_parts * kWordBits * kWordBits);
        }

        char *inverse_bit_matrix = new char[m_num_data_parts * m_num_data_parts
                                            * kWordBits * kWordBits];
        _InvertBitMatrix(tmp_bit_matrix, inverse_bit_matrix,
                        m_num_data_parts * kWordBits);

        iter = decoding_bit_matrix;
        for (int i = 0; i < num_erased_data_parts; i++) {
            memcpy(iter, inverse_bit_matrix + m_num_data_parts * kWordBits
                    * kWordBits * rowid_to_partidx[m_num_data_parts + i],
                    sizeof(char) * m_num_data_parts * kWordBits * kWordBits); // NOLINT
            iter += m_num_data_parts*kWordBits*kWordBits;
        }

        delete[] tmp_bit_matrix;
        delete[] inverse_bit_matrix;
    }

    /* Next, here comes the hard part.  For each coding node that needs
       to be decoded, you start by putting its rows of the distribution
       matrix into the decoding matrix.  If there were no failed data
       nodes, then you're done.  However, if there have been failed
       data nodes, then you need to modify the columns that correspond
       to the data nodes.  You do that by first zeroing them.  Then
       whereever there is a one in the distribution matrix, you XOR
       in the corresponding row from the failed data node's entry in
       the decoding matrix.  The whole process kind of makes my head
       spin, but it works.
     */

    for (int k = 0; k < num_erased_code_parts; k++) {
        int code_part_idx = rowid_to_partidx[k + num_erased_data_parts + m_num_data_parts]
                                - m_num_data_parts;

        char *iter = decoding_bit_matrix + m_num_data_parts * kWordBits * kWordBits
                    * (num_erased_data_parts + k);
        memcpy(iter,
                m_encoding_bit_matrix + code_part_idx * m_num_data_parts * kWordBits * kWordBits,
                sizeof(char) * m_num_data_parts * kWordBits * kWordBits); // NOLINT

        for (int i = 0; i < m_num_data_parts; i++) {
            if (rowid_to_partidx[i] != i) {
                for (int j = 0; j < kWordBits; j++) {
                    bzero(iter + j * m_num_data_parts * kWordBits + i * kWordBits,
                            sizeof(char) * kWordBits); // NOLINT
                }
            }
        }

        /* There's the yucky part */
        int index = code_part_idx * m_num_data_parts * kWordBits * kWordBits;
        for (int i = 0; i < m_num_data_parts; i++) {
            if (rowid_to_partidx[i] != i) {
                char *p1 = decoding_bit_matrix + (partidx_to_rowid[i] - m_num_data_parts)
                            *m_num_data_parts*kWordBits*kWordBits;
                for (int j = 0; j < kWordBits; j++) {
                    char *p2 = iter + j*m_num_data_parts*kWordBits;
                    for (int m = 0; m < kWordBits; m++) {
                        if (m_encoding_bit_matrix[index
                                                    + j * m_num_data_parts * kWordBits
                                                    + i * kWordBits + m]) {
                            for (int n = 0; n < m_num_data_parts * kWordBits; n++) {
                                p2[n] = p2[n] ^ p1[n + m * m_num_data_parts * kWordBits];
                            }
                        }
                    }
                }
            }
        }
    }

    // Generate decoding schedule
    int **decoding_schedule = _BitMatrixToSchedule(m_num_data_parts,
                                                    num_erased_data_parts + num_erased_code_parts,
                                                    decoding_bit_matrix);

    // do decoding
    _DoScheduleOperations(decoding_schedule, ptrs, size);

    delete[] decoding_bit_matrix;
    _FreeSchedule(decoding_schedule);
}

}
}
