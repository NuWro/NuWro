#pragma once

#include <cstddef>
#include <fcntl.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

template <typename T> class mmapio {
private:
  T *data;
  size_t m_size;
  int fd;

public:
  mmapio(const std::string filename, bool allow_write, size_t count = 0)
      : mmapio(filename.c_str(), allow_write, count) {}
  mmapio(const char *filename, bool allow_write, size_t count = 0) {
    try {
      auto size = count * sizeof(T);
      auto openmode = allow_write ? O_RDWR : O_RDONLY;
      if (allow_write)
        openmode |= O_CREAT;
      fd = open(filename, openmode, 0644);
      if (fd == -1) {
        throw std::runtime_error("open failed");
      }
      struct stat64 file_stat;
      if (fstat64(fd, &file_stat) == -1) {
        throw std::runtime_error("fstat failed");
      }
      if (size) {
        if (file_stat.st_size < size) {
          if (allow_write) {
            if (ftruncate(fd, size) == -1) {
              throw std::runtime_error("ftruncate failed");
            }
          } else {
            throw std::runtime_error("file size is smaller than requested, "
                                     "while open in read-only mode");
          }
        }
      } else {
        size = file_stat.st_size;
      }
      if (size == 0) {
        throw std::runtime_error("size is 0");
      }
      m_size = size;
      auto page_mode = allow_write ? (PROT_READ | PROT_WRITE) : PROT_READ;
      data =
          static_cast<T *>(mmap(nullptr, size, page_mode, MAP_SHARED, fd, 0));
      if (data == MAP_FAILED) {
        throw std::runtime_error("mmap failed");
      }
    } catch (std::runtime_error &e) {
      if (data != MAP_FAILED) {
        munmap(data, m_size);
      }
      if (fd != -1) {
        close(fd);
      }
      throw e;
    }
  }
  mmapio(const mmapio &) = delete;
  mmapio(mmapio &&other) {
    data = other.data;
    m_size = other.m_size;
    fd = other.fd;
    other.data = nullptr;
    other.m_size = 0;
    other.fd = -1;
  }
  ~mmapio() {
    if (fd == -1) { // in case this is a moved object
      return;
    }
    if (munmap(data, m_size) == -1) {
      std::cerr << "munmap failed" << std::endl
                << "errno: " << errno << std::endl;
    }
    if (close(fd) == -1) {
      std::cerr << "close failed" << std::endl;
    }
  }
  T &operator[](size_t index) const { return data[index]; }
  size_t count() const { return m_size / sizeof(T); }
  T *begin() const { return data; }
  T *end() const { return data + count(); }
};