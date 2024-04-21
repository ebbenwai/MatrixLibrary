#include "s21_matrix.h"

void s21_copy_matrix(matrix_t *A, matrix_t *B) {
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      B->matrix[i][j] = A->matrix[i][j];
    }
  }
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (result == NULL || rows <= 0 || columns <= 0) return INC_MX;

  int code = OK;

  if ((result->matrix = (double **)calloc(rows, sizeof(double *)))) {
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
      if (result->matrix[i] == NULL) code = INC_MX;
    }
  } else
    code = INC_MX;

  if (code == OK) {
    result->rows = rows;
    result->columns = columns;
  }

  return code;
}

void s21_remove_matrix(matrix_t *A) {
  if (A == NULL || A->rows < 0 || A->columns < 0) {
    return;
  }
  for (int i = 0; i < A->rows; i++) {
    if (A->matrix[i] != NULL) {
      free(A->matrix[i]);
    }
  }
  if (A->matrix != NULL) {
    free(A->matrix);
  }
  A->matrix = NULL;
  A->rows = 0;
  A->columns = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int code = SUCCESS;
  if (A == NULL || B == NULL || (A->rows != B->rows) ||
      (A->columns != B->columns))
    code = FAILURE;

  for (int i = 0; i < A->rows && code == SUCCESS; i++) {
    for (int j = 0; j < A->columns && code == SUCCESS; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
        code = FAILURE;
      }
    }
  }

  return code;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return INC_MX;

  int flag = OK;

  if (A->rows != B->rows || A->columns != B->columns) {
    flag = CALC_ERR;
  } else {
    flag = s21_create_matrix(A->rows, A->columns, result);

    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        if (isfinite(result->matrix[i][j])) {
          result->matrix[i][j] = result->matrix[i][j];
        } else {
          flag = CALC_ERR;
        }
      }
    }
  }

  return flag;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return INC_MX;

  int flag = OK;

  if (A->rows != B->rows || A->columns != B->columns) {
    flag = CALC_ERR;
  } else {
    flag = s21_create_matrix(A->rows, A->columns, result);

    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        if (isfinite(result->matrix[i][j])) {
          result->matrix[i][j] = result->matrix[i][j];
        } else {
          flag = CALC_ERR;
        }
      }
    }
  }

  return flag;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return INC_MX;

  int code = OK;

  if (A->columns != B->rows)
    code = CALC_ERR;
  else
    code = s21_create_matrix(A->rows, B->columns, result);

  for (int i = 0; i < A->rows && code == OK; i++) {
    for (int k = 0; k < B->columns; k++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][k] += A->matrix[i][j] * B->matrix[j][k];
      }
    }
  }

  return code;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (A == NULL || result == NULL) return INC_MX;
  int flag = OK;
  flag = s21_create_matrix(A->rows, A->columns, result);

  for (int i = 0; i < A->rows && flag == OK; ++i) {
    for (int j = 0; j < A->columns && flag == OK; ++j) {
      result->matrix[i][j] = number * A->matrix[i][j];
      if (isfinite(result->matrix[i][j])) {
        result->matrix[i][j] = result->matrix[i][j];
      } else {
        flag = CALC_ERR;
      }
    }
  }

  return flag;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int flag = OK;

  if (A == NULL || A->columns <= 0 || A->rows <= 0 || result == NULL)
    return INC_MX;
  s21_create_matrix(A->columns, A->rows, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[j][i] = A->matrix[i][j];
    }
  }

  return flag;
}

int s21_determinant(matrix_t *A, double *result) {
  int flag = OK;
  if (A == NULL || result == NULL || A->columns <= 0 || A->rows <= 0)
    return INC_MX;
  if (A->rows != A->columns) return CALC_ERR;

  if (A->rows == 1) {
    *result = A->matrix[0][0];
  } else if (A->rows == 2) {
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[1][0] * A->matrix[0][1];
  } else {
    double det = 1.0;

    for (int k = 0; k < A->rows; k++) {
      int max_row = k;
      for (int i = k + 1; i < A->columns; i++) {
        if (fabs(A->matrix[i][k]) > fabs(A->matrix[max_row][k])) {
          max_row = i;
        }
      }

      if (max_row != k) {
        for (int j = 0; j < A->rows; j++) {
          double temp = A->matrix[k][j];
          A->matrix[k][j] = A->matrix[max_row][j];
          A->matrix[max_row][j] = temp;
        }
        det *= -1;
      }

      det *= A->matrix[k][k];

      for (int i = k + 1; i < A->rows; i++) {
        double ratio = A->matrix[i][k] / A->matrix[k][k];
        for (int j = k; j < A->columns; j++) {
          A->matrix[i][j] -= ratio * A->matrix[k][j];
        }
      }
    }

    *result = det;
  }

  return flag;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL || A->columns <= 0 || A->rows <= 0)
    return INC_MX;

  if (A->rows != A->columns || A->rows <= 0 || A->columns <= 0) return CALC_ERR;

  int flag = OK;
  flag = s21_create_matrix(A->rows, A->columns, result);

  matrix_t minor;
  s21_create_matrix(A->rows - 1, A->columns - 1, &minor);

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      int mi = 0;
      for (int m = 0; m < A->rows; m++) {
        if (m == i) {
          continue;
        }

        int mj = 0;
        for (int n = 0; n < A->columns; n++) {
          if (n == j) {
            continue;
          }
          minor.matrix[mi][mj] = A->matrix[m][n];

          mj++;
        }
        mi++;
      }

      double minor_det = 0;
      s21_determinant(&minor, &minor_det);

      result->matrix[i][j] = minor_det;

      result->matrix[i][j] *= pow(-1, (i + j));
    }
  }

  s21_remove_matrix(&minor);
  return flag;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int flag = OK;

  if (A == NULL || result == NULL || A->columns == 0 || A->rows == 0)
    return INC_MX;

  matrix_t B;
  s21_create_matrix(A->rows, A->columns, &B);

  s21_copy_matrix(A, &B);

  double determinate = 0;

  s21_determinant(&B, &determinate);

  if (determinate != 0 && A->rows == 1 && A->columns == 1 && flag == OK) {
    flag = s21_create_matrix(A->columns, A->rows, result);
    result->matrix[0][0] = 1 / A->matrix[0][0];
  } else if (determinate != 0 && flag == OK) {
    double det = 1.0 / determinate;

    matrix_t calc_complements;
    matrix_t transpose;

    flag = s21_calc_complements(A, &calc_complements);
    flag = s21_transpose(&calc_complements, &transpose);
    flag = s21_mult_number(&transpose, det, result);

    s21_remove_matrix(&calc_complements);
    s21_remove_matrix(&transpose);

  } else {
    flag = CALC_ERR;
  }
  s21_remove_matrix(&B);

  return flag;
}
