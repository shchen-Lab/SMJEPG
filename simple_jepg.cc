// Copyright (c) 2020 smarsufan. All Rights Reserved.

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <cstring>
#include <cassert>
#include <cmath>
#include "simple_jepg.h"

static uint8_t zip_table_[64] = {
  0,  1,  5,  6,  14, 15, 27, 28,
  2,  4,  7,  13, 16, 26, 29, 42,
  3,  8,  12, 17, 25, 30, 41, 43,
  9,  11, 18, 24, 31, 40, 44, 53,
  10, 19, 23, 32, 39, 45, 52, 54,
  20, 22, 33, 38, 46, 41, 55, 60,
  21, 34, 37, 47, 50, 56, 59, 61,
  35, 36, 48, 49, 57, 58, 62, 63,
};

static float MtxIDCT[64] = {
  0.3536,    0.4904,    0.4619,    0.4157,    0.3536,    0.2778,    0.1913,    0.0975,
  0.3536,    0.4157,    0.1913,   -0.0975,   -0.3536,   -0.4904,   -0.4619,   -0.2778,
  0.3536,    0.2778,   -0.1913,   -0.4904,   -0.3536,    0.0975,    0.4619,    0.4157,
  0.3536,    0.0975,   -0.4619,   -0.2778,    0.3536,    0.4157,   -0.1913,   -0.4904,
  0.3536,   -0.0975,   -0.4619,    0.2778,    0.3536,   -0.4157,   -0.1913,    0.4904,
  0.3536,   -0.2778,   -0.1913,    0.4904,   -0.3536,   -0.0975,    0.4619,   -0.4157,
  0.3536,   -0.4157,    0.1913,    0.0975,   -0.3536,    0.4904,   -0.4619,    0.2778,
  0.3536,   -0.4904,    0.4619,   -0.4157,    0.3536,   -0.2778,    0.1913,   -0.0975,
};

static float MtxDCT[64] = {
  0.3536,    0.3536,    0.3536,    0.3536,    0.3536,    0.3536,    0.3536,    0.3536,
  0.4904,    0.4157,    0.2778,    0.0975,   -0.0975,   -0.2778,   -0.4157,   -0.4904,
  0.4619,    0.1913,   -0.1913,   -0.4619,   -0.4619,   -0.1913,    0.1913,    0.4619,
  0.4157,   -0.0975,   -0.4904,   -0.2778,    0.2778,    0.4904,    0.0975,   -0.4157,
  0.3536,   -0.3536,   -0.3536,    0.3536,    0.3536,   -0.3536,   -0.3536,    0.3536,
  0.2778,   -0.4904,    0.0975,    0.4157,   -0.4157,   -0.0975,    0.4904,   -0.2778,
  0.1913,   -0.4619,    0.4619,   -0.1913,   -0.1913,    0.4619,   -0.4619,    0.1913,
  0.0975,   -0.2778,    0.4157,   -0.4904,    0.4904,   -0.4157,    0.2778,   -0.0975
};

template<typename T, typename... Ts>
std::unique_ptr<T> make_unique(Ts&&... params) {
    return std::unique_ptr<T>(new T(std::forward<Ts>(params)...));
}

int64_t timeus() {
  struct timeval tv;  
  gettimeofday(&tv,NULL);
  return tv.tv_sec * 1000000 + tv.tv_usec; 
}

std::unique_ptr<std::vector<uint8_t>> Read(const std::string &path) {
  std::ifstream fb(path);
  if (!fb.is_open()) {
    return nullptr;
  }

  fb.seekg(0, fb.end);
  size_t size = fb.tellg();
  fb.seekg(0, fb.beg);

  std::unique_ptr<std::vector<uint8_t>> vec = make_unique<std::vector<uint8_t>>(size);
  fb.read(reinterpret_cast<char *>(vec->data()), size);
  fb.close();

  return std::move(vec);
}

template <typename T, typename V, typename R>
void matmul(T *A, V *B, R *C) {
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      R sum = 0;
      for (int k = 0; k < 8; ++k) {
        sum += A[i * 8 + k] * B[k * 8 + j];
      }
      C[i * 8 + j] = round(sum);
    }
  }
}

JPEGDecoder::JPEGDecoder(const std::string path) {
  raw_ = Read(path);
  raw_data_ = raw_->data();
  Decode();
}

int pack(uint8_t *ptr, int size, bool isLittleEnd) {
  // TODO: Decode little endian and big endian.
  return -1;
}

void JPEGDecoder::Decode() {
  GetQuantTable();
  GetSize();
  GetHuffmanTable();

  DecodeData();
}

void JPEGDecoder::GetQuantTable() {
  quant_tables_[0] = raw_data_ + 25;
  quant_tables_[1] = raw_data_ + 25 + 64 + 5;

#ifndef NDEBUG
  std::cout << "quant Y: " << std::endl;
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      std::cout << int(quant_tables_[0][i * 8 + j]) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "quant C: " << std::endl;
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      std::cout << int(quant_tables_[1][i * 8 + j]) << " ";
    }
    std::cout << std::endl;
  }
#endif
}

void JPEGDecoder::GetSize() {
  uint8_t *height_ptr = raw_data_ + 163;
  uint8_t *width_ptr = raw_data_ + 165;
  height_ = (height_ptr[0] << 8) | height_ptr[1];
  width_ = (width_ptr[0] << 8) | width_ptr[1];

#ifndef NDEBUG
  std::cout << "height ... " << height_ << " width ... " << width_ << std::endl;
#endif
}

void JPEGDecoder::GetHuffmanTable() {
  uint8_t *ptr = raw_data_ + 181;
  ptr = DecodeHuffmanTable(ptr);

  ptr += 4;
  ptr = DecodeHuffmanTable(ptr);

  ptr += 4;
  ptr = DecodeHuffmanTable(ptr);

  ptr += 4;
  ptr = DecodeHuffmanTable(ptr);

#ifndef NDEBUG
  huffman.show();
#endif
}

uint8_t *JPEGDecoder::DecodeHuffmanTable(uint8_t *ptr) {
  uint8_t table_idx = ptr[0];
  auto &table = huffman[table_idx];

  uint8_t *num_keys_ptr = ptr + 1;
  uint8_t *val_ptr = num_keys_ptr + 16;

  uint16_t pre_bits_value = 0;
  for (int idx = 0; idx < 16; ++idx) {
    pre_bits_value <<= 1;

    int high_key = (idx + 1) << 16;

    uint8_t num = num_keys_ptr[idx];
    for (uint8_t n = 0; n < num; ++n) {
      uint8_t val = (val_ptr++)[0];

      int low_key = (n + pre_bits_value);
      table[high_key | low_key] = val;
    }

    pre_bits_value += num;
  }

  return val_ptr;
}

void JPEGDecoder::DecodeData() {
  im_data_ = raw_data_ + 623;
  bit_value_ = im_data_[0];

  int h_block = (height_ + 16 - 1) / 16;
  int w_block = (width_ + 16 - 1) / 16;

  float matrix[6][64];

  data_ = new uint8_t[height_ * width_ * 3];

  for (int h_block_id = 0; h_block_id < h_block; ++h_block_id) {
    for (int w_block_id = 0; w_block_id < w_block; ++w_block_id) {
      memset(matrix, 0, sizeof(matrix));

      // Decode matrix Y
      DecodeMatrix(matrix[0], 0, quant_tables_[0], dc[0]);
      DecodeMatrix(matrix[1], 0, quant_tables_[0], dc[0]);
      DecodeMatrix(matrix[2], 0, quant_tables_[0], dc[0]);
      DecodeMatrix(matrix[3], 0, quant_tables_[0], dc[0]);

      // Decode matrix Cb
      DecodeMatrix(matrix[4], 1, quant_tables_[1], dc[1]);
      // Decode matrix Cr
      DecodeMatrix(matrix[5], 1, quant_tables_[1], dc[2]);

      // YUV2RGB
      {
        int h_offset = h_block_id * 16;
        int w_offset = w_block_id * 16;

        YUV2RGB(data_, matrix[0], matrix[4], matrix[5], std::min(8, height_ - h_offset), std::min(8, width_ - w_offset), h_offset, w_offset);
        YUV2RGB(data_, matrix[1], matrix[4] + 4, matrix[5] + 4, std::min(8, height_ - h_offset), std::min(8, width_ - w_offset + 8), h_offset, w_offset + 8);
        YUV2RGB(data_, matrix[2], matrix[4] + 4 * 8, matrix[5] + 4 * 8, std::min(8, height_ - h_offset - 8), std::min(8, width_ - w_offset), h_offset + 8, w_offset);
        YUV2RGB(data_, matrix[3], matrix[4] + 4 * 8 + 4, matrix[5] + 4 * 8 + 4, std::min(8, height_ - h_offset - 8), std::min(8, width_ - w_offset + 8), h_offset + 8, w_offset + 8);
      }
    }
  }
}

void JPEGDecoder::DecodeMatrix(float *matrix, uint8_t table_idx, uint8_t *quant, int &dc) {
  auto &table = huffman[table_idx];
  uint8_t len = FindInHuffmanTable(table);
  int value = GetValue(len);
  dc += value;
  matrix[0] = dc;

  auto &ac_table = huffman[table_idx | 16];
  for (int idx = 1; idx < 64; ++idx) {
    uint8_t len = FindInHuffmanTable(ac_table);
    if (len == 0) break;

    value = GetValue(len & 0xf);
    idx += (len >> 4);
    matrix[idx] = value;
  }

  for (int idx = 0; idx < 64; ++idx) {
    matrix[idx] *= quant[idx];
  }

  unzip(matrix);
  IDCT(matrix);
}

uint8_t JPEGDecoder::FindInHuffmanTable(std::unordered_map<int, uint8_t> &table) {
  int val = 0;
  for (int idx = 0; idx < 16; ++idx) {
    val <<= 1;
    val += NextBit();
    int r_val = ((idx + 1) << 16) | val;
    if (table.find(r_val) != table.end()) {
      return table[r_val];
    }
  }
  assert(0);
}

int JPEGDecoder::NextBit() {
  if (bit_mask_ == 0x00) {
    bit_mask_ = 0x80;
    // bit_offset = 0;

    bit_value_ = (++im_data_)[0];
    if (bit_value_ == 0xFF) {
      bit_value_ = (++im_data_)[0];

      switch (bit_value_) {
        case 0x00:
          bit_value_ = 0xFF;
          break;

        case 0xD7:
          dc[0] = dc[1] = dc[2] = 0;
          break;

        default:
          break;
      }
    }
  }

  int val = bit_mask_ & bit_value_;
  bit_mask_ >>= 1;
  // bit_offset += 1;
  return val > 0 ? 1 : 0;
}

int JPEGDecoder::GetValue(uint8_t len) {
  // if (len == 0) return 0;

  // int tag = -2 * NextBit() + 1;

  // int value = 0;
  // for (int idx = 1; idx < len; ++idx) {
  //   value = (value << 1) | NextBit();
  // }
  // value *= tag;
  // return value;
    int value = 0;
    for (int idx = 0; idx < len; ++idx) {
      value = (value << 1) + NextBit();
    }
    return value >= pow(2, len - 1) ? value : value - pow(2, len) + 1;
}

template <typename T>
void JPEGDecoder::unzip(T *matrix) {
  memcpy(unzip_matrix_, matrix, sizeof(T) * 64);
  for (int idx = 0; idx < 64; ++idx) {
    matrix[idx] = unzip_matrix_[zip_table_[idx]];
  }
}

template <typename T>
void JPEGDecoder::IDCT(T *matrix) {
  matmul(MtxIDCT, matrix, unzip_matrix_);
  matmul(unzip_matrix_, MtxDCT, matrix);
}

template <typename T>
void JPEGDecoder::YUV2RGB(uint8_t *data, T *Y, T *Cr, T *Cb, int h, int w, int rh, int rw) {
  for (int i = 0; i < h; ++i) {
    for (int j = 0; j < w; ++j) {
      uint8_t *rgb = data + ((i + rh) * width_ + j + rw) * 3;

      T y = Y[i * 8 + j];
      T cr = Cr[i / 2 * 8 + j / 2];
      T cb = Cb[i / 2 * 8 + j / 2];

      T R = y + 1.402f * cr + 128.f;
      T G = y - 0.34414f * cb - 0.71414f * cr + 128.f;
      T B = y + 1.772f * cb + 128.f;

      R = std::min(std::max(R, 0.f), 255.f);
      G = std::min(std::max(G, 0.f), 255.f);
      B = std::min(std::max(B, 0.f), 255.f);

      rgb[0] = R;
      rgb[1] = G;
      rgb[2] = B;
    }
  }
}
