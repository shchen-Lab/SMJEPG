// Copyright (c) 2020 smarsufan. All Rights Reserved.

#pragma once
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <iostream>

int64_t timeus();

typedef unsigned const char *BufPtr;

static std::unique_ptr<std::vector<uint8_t>> Read(const std::string &path);

class HuffmanTable {
 public:
  std::unordered_map<int, uint8_t> &operator[](uint8_t idx) {
    // [0, 1, 16, 17] -> [0, 1, 2, 3]
    // std::cout << "Get huffman Table " << ((idx >> 3) | (idx & 0x01)) << std::endl;
    return tables[(idx >> 3) | (idx & 0x01)];
  }

  void show() {
    std::cout << "Y_DC: " << std::endl;
    for (auto iter = tables[0].begin(); iter != tables[0].end(); ++iter) {
      std::cout << iter->first << ": " << int(iter->second) << std::endl;
    }

    std::cout << "Y_AC: " << std::endl;
    for (auto iter = tables[1].begin(); iter != tables[1].end(); ++iter) {
      std::cout << iter->first << ": " << int(iter->second) << std::endl;
    }

    std::cout << "C_DC: " << std::endl;
    for (auto iter = tables[2].begin(); iter != tables[2].end(); ++iter) {
      std::cout << iter->first << ": " << int(iter->second) << std::endl;
    }

    std::cout << "C_AC: " << std::endl;
    for (auto iter = tables[3].begin(); iter != tables[3].end(); ++iter) {
      std::cout << iter->first << ": " << int(iter->second) << std::endl;
    }
  }
  
 private:
  std::unordered_map<int, uint8_t> tables[4];
};

class JPEGDecoder {
 public:
  explicit JPEGDecoder(const std::string path);

  int h() { return height_; }
  
  int w() { return width_; }

  uint8_t *data() { return data_; }

 private:
  int pack(uint8_t *ptr, int size, bool isLittleEnd);
 
  void Decode();

  void GetQuantTable();

  void GetSize();
  // Get four huffman table
  void GetHuffmanTable();

  // decode one huffman table
  uint8_t *DecodeHuffmanTable(uint8_t *ptr);

  void DecodeData();

  void DecodeMatrix(float *matrix, uint8_t table_idx, uint8_t *quant, int &dc);

  int NextBit();

  uint8_t FindInHuffmanTable(std::unordered_map<int, uint8_t> &table);

  int GetValue(uint8_t len);

  template <typename T>
  void unzip(T *matrix);

  template <typename T>
  void IDCT(T *matrix);

  template <typename T>
  void YUV2RGB(uint8_t *data, T *Y, T *Cr, T *Cb, int h, int w, int rh, int rw);
  
  bool findJPEGheader(BufPtr *start, uint32_t *len, uint8_t marker);
  bool decodeJPEGfile(BufPtr *start, uint32_t *len, BufPtr *qtable0, BufPtr *qtable1);
  void skipScanBytes(BufPtr *start);
  void nextJpegBlock(BufPtr *bytes);
 private:
  std::unique_ptr<std::vector<uint8_t>> raw_;

  HuffmanTable huffman;

  // {quant table Y, quant table C}
  uint8_t *quant_tables_[2] = {nullptr, nullptr};

  int dc[3] = {0, 0, 0};

  float unzip_matrix_[64];

  uint8_t  *raw_data_;

  uint8_t *im_data_;

  uint8_t *data_{nullptr};

  int height_{0};

  int width_{0};

  uint8_t bit_mask_{0x80};

  // uint8_t bit_offset{0};

  uint8_t bit_value_{0};

  bool is_end_{false};
};
