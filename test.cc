// Copyright (c) 2020 smarsufan. All Rights Reserved.

#include <iostream>
#include <fstream>
#include "simple_jepg.h"

int main(int argv, char *args[]) {
  if (argv != 2) {
    return -1;
  }

  std::string path(args[1]);
  auto t1 = timeus();
  JPEGDecoder decoder(path);
  auto t2 = timeus();

  std::cout << "time ... " << (t2 - t1) / 1000. << " ms" << std::endl;

  std::ofstream fb("bird.txt");
  for (int idx = 0; idx < decoder.h() * decoder.w() * 3; ++idx) {
    fb << int(decoder.data()[idx]) << std::endl;
  }

  return 0;
}
