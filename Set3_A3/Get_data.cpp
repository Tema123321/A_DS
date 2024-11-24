// ConsoleApplication11.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <algorithm>

class Random {
 private:
  int64_t value;
  int64_t A = 48721;
  int64_t M = 2147483647;

 public:
  Random() { value = 1; }
  Random(int64_t seed) { value = seed; }
  int64_t get_int_random(int64_t start, int64_t end) {
    value = (A * value) % M;
    return value % (end - start + 1) + start;
  }
};

void insertionSort(std::vector<int64_t> &A, int64_t start, int64_t end) {
  for (int64_t i = start; i <= end; ++i) {
    int64_t ptr = i;
    int64_t x = A[i];
    while (ptr > start && A[ptr - 1] > x) {
      A[ptr] = A[ptr - 1];
      A[ptr-- - 1] = x;
    }
  }
}

int64_t partition(std::vector<int64_t> &A, int64_t start, int64_t end, int64_t pivot_id) {
  int64_t pivot = A[pivot_id];
  std::swap(A[start], A[pivot_id]);
  int64_t left_ptr = start + 1;
  int64_t right_ptr = end;
  while (true) {
    while (left_ptr <= end && A[left_ptr] <= pivot) {
      ++left_ptr;
    }
    while (right_ptr >= start && A[right_ptr] > pivot) {
      --right_ptr;
    }
    if (left_ptr < right_ptr) {
      std::swap(A[left_ptr], A[right_ptr]);
    } else {
      break;
    }
  }
  std::swap(A[start], A[right_ptr]);
  return right_ptr;
}

void heapify(std::vector<int64_t> &A, int64_t now, int64_t start, int64_t end) {
 while (2 * (now - start) + start + 1 <= end &&
          A[now] < A[2 * (now - start) + start + 1] ||
      2 * (now - start) + start + 2 <= end &&
          A[now] < A[2 * (now - start) + start + 2]) {
    if (2 * (now - start) + start + 2 <= end &&
        A[2 * (now - start) + start + 2] > A[2 * (now - start) + start + 1]) {
      std::swap(A[now], A[2 * (now - start) + start + 2]);
      now = 2 * (now - start) + start + 2;
      continue;
    }
    std::swap(A[now], A[2 * (now - start) + start + 1]);
    now = 2 * (now - start) + start + 1;
  }
}

void buildmaxheap(std::vector<int64_t> &A, int64_t start, int64_t end) {
  for (int64_t i = (end - start - 1) / 2 + start; i >= start; --i) {
    heapify(A, i, start, end);
  }
}

void heapsort(std::vector<int64_t> &A, int64_t start, int64_t end) {
  buildmaxheap(A, start, end);
  for (int64_t i = end; i > start; --i) {
    std::swap(A[start], A[i]);
    heapify(A, start, start, i - 1);
  }
}

void quicksort(std::vector<int64_t> &A, int64_t start, int64_t end, int64_t deep, Random &r) {
  if (end - start + 1 <= 15) {
    return insertionSort(A, start, end);
  }
  if (std::pow(2, deep) >= A.size() * A.size()) {
    return heapsort(A, start, end);
  }
  int64_t pivot_id = r.get_int_random(start, end);
  pivot_id = partition(A, start, end, pivot_id);
  quicksort(A, start, pivot_id - 1, deep + 1, r);
  quicksort(A, pivot_id + 1, end, deep + 1, r);
}

void quicksortN(std::vector<int64_t> &A, int64_t start, int64_t end, Random &r) {
  if (end - start + 1 <= 1) {
    return;
  }
  int64_t pivot_id = r.get_int_random(start, end);
  pivot_id = partition(A, start, end, pivot_id);
  quicksortN(A, start, pivot_id - 1, r);
  quicksortN(A, pivot_id + 1, end, r);
}

void introSort(std::vector<int64_t> &A) {
  Random r(2645);
  quicksort(A, 0, A.size() - 1, 1, r);
}

void quickSort(std::vector<int64_t> &A) {
  Random r(2645);
  quicksortN(A, 0, A.size() - 1, r);
}

class ArrayGenerator {
 private:
  Random r;

 public:
  ArrayGenerator(int64_t seed) { r = Random(seed); }

  std::vector<int64_t> RandomVector(size_t size, int64_t start, int64_t end) {
    std::vector<int64_t> A(size);
    for (size_t i = 0; i < size; ++i) {
      A[i] = r.get_int_random(start, end);
    }
    return A;
  }

  std::vector<int64_t> AntiSortVector(size_t size, int64_t start, int64_t end) {
    std::vector<int64_t> A = RandomVector(size, -end, -start);
    std::sort(A.begin(), A.end());
    for (size_t i = 0; i < size; ++i) {
      A[i] = -A[i];
    }
    return A;
  }

  std::vector<int64_t> ShuffleSortVector(size_t size, int64_t cnt_shuffles,
                                         int64_t start, int64_t end) {
    std::vector<int64_t> A = RandomVector(size, start, end);
    std::sort(A.begin(), A.end());
    for (int64_t i = 0; i < cnt_shuffles; ++i) {
      int64_t x = r.get_int_random(0, size - 1);
      int64_t y = r.get_int_random(0, size - 1);
      std::swap(A[x], A[y]);
    }
    return A;
  }
};

class SortTester {
 public:
  double TestSort(const std::vector<int64_t> &A, int64_t l, int64_t r) {
    std::vector<int64_t> B(r - l + 1);
    for (int64_t i = 0; i <= r - l; ++i) {
      B[i] = A[l + i];
    }
    auto start = std::chrono::high_resolution_clock::now();
    std::sort(B.begin(), B.end());
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    return duration.count();
  }

  double TestQuickSort(const std::vector<int64_t> &A, int64_t l, int64_t r, Random rand) {
    std::vector<int64_t> B(r - l + 1);
    for (int64_t i = 0; i <= r - l; ++i) {
      B[i] = A[l + i];
    }
    auto start = std::chrono::high_resolution_clock::now();
    quickSort(B);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    return duration.count();
  }

  double TestQuickInsertHeapSort(const std::vector<int64_t> &A, int64_t l, int64_t r, Random rand) {
    std::vector<int64_t> B(r - l + 1);
    for (int64_t i = 0; i <= r - l; ++i) {
      B[i] = A[l + i];
    }
    auto start = std::chrono::high_resolution_clock::now();
    introSort(B);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    return duration.count();
  }
};

void print(const std::vector<int64_t> &v) {
  for (size_t i = 0; i < v.size(); ++i) {
    std::cout << v[i] << " ";
  }
  std::cout << "\n";
}

int main() {
  ArrayGenerator generator(67798);
  SortTester st;
  std::string file = "C:\\Users\\kas11\\Downloads\\A3_data.txt";
  std::ofstream output(file, std::ios::out);

  int64_t n = 1;
  std::vector<std::vector<double>> data(501);
  for (int64_t i = 1; i < 501; ++i) {
    data[i] = std::vector<double>(9);
  }

  Random r(143553);
  std::vector<Random> randoms(n);
  for (int64_t i = 0; i < n; ++i) {
    randoms[i] = Random(r.get_int_random(1, 1000000));
  }

  for (int64_t iter = 0; iter < n; ++iter) {
    std::vector<int64_t> A_rand = generator.RandomVector(50000, 0, 45000);
    std::vector<int64_t> A_anti_sort =
        generator.AntiSortVector(50000, 0, 45000);
    std::vector<int64_t> A_shuffle =
        generator.ShuffleSortVector(50000, 5000, 0, 45000);
    for (int64_t size = 100; size <= 50000; size += 100) {
      data[size / 100][0] +=
          st.TestQuickSort(A_rand, 0, size - 1, randoms[iter]);
      data[size / 100][1] +=
          st.TestQuickSort(A_anti_sort, 0, size - 1, randoms[iter]);
      data[size / 100][2] +=
          st.TestQuickSort(A_shuffle, 0, size - 1, randoms[iter]);
      data[size / 100][3] +=
          st.TestQuickInsertHeapSort(A_rand, 0, size - 1, randoms[iter]);
      data[size / 100][4] +=
          st.TestQuickInsertHeapSort(A_anti_sort, 0, size - 1, randoms[iter]);
      data[size / 100][5] +=
          st.TestQuickInsertHeapSort(A_shuffle, 0, size - 1, randoms[iter]);
      data[size / 100][6] += st.TestSort(A_rand, 0, size - 1);
      data[size / 100][7] += st.TestSort(A_anti_sort, 0, size - 1);
      data[size / 100][8] += st.TestSort(A_shuffle, 0, size - 1);
    }
  }

  for (int64_t i = 1; i < 501; ++i) {
    for (int64_t j = 0; j < 9; ++j) {
      data[i][j] /= n;
    }
  }

  for (int64_t i = 1; i < 501; ++i) {
    output << i * 100 << ";";
    for (int64_t j = 0; j < 9; ++j) {
      output << data[i][j] << ";";
    }
    output << std::endl;
  }
  output.close();
}