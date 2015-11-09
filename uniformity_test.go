package bloom

import (
	"testing"
	"github.com/spaolacci/murmur3"
	"fmt"
	"math"
	"hash/fnv"
)

// I set this up as a type so that I can pass functions into a tester.
type locations func(data []byte, k, m uint32) []uint

// This mimics the implementation of f.Add(), but without the filter.
func bloomLocations(data []byte, k, m uint32) []uint {
	locs := make([]uint, k)

	h := baseHashes(data)

	// Methods like Add() range from ii = 0 to k
	for i := uint32(0); i < k; i++ {
		// This exactly re-implements location()
		ii := uint64(i)
		locs[i] = uint((h[ii%2] + ii*h[2+(((ii+(ii%2))%4)/2)]) % uint64(m))
	}

	return locs
}

// This does the same, but prior to changeset c3a7944627e42592b577b8014c6bed5351100825
func wfnvLocations(data []byte, k, m uint32) []uint {
	locs := make([]uint, k)

	h := baseHashes(data)

	// Methods like Add() range from ii = 0 to k
	for i := uint32(0); i < k; i++ {
		// This exactly re-implements location()
		ii := uint64(i)
		locs[i] = uint((h[0] + ii*h[1]) % uint64(m))
	}

	return locs
}

// This uses the native FNV-1 implementation.
func fnvLocations(data []byte, k, m uint32) []uint {
	locs := make([]uint, k)

	hash := fnv.New64()
	hash.Write(data)

	// Sadly, because this is a 64 bit hasher, we can't split
	// it and keep 64-bit addressing.
	h0 := uint64(uint32(hash.Sum64()))
	h1 := uint64(uint32(hash.Sum64() >> 32))

	for i := uint32(0); i < k; i++ {
		ii := uint64(i)
		locs[i] = uint((h0 + ii*h1) % uint64(m))
	}

	return locs
}

// Same thing, but with FNV-1a
func fnvaLocations(data []byte, k, m uint32) []uint {
	locs := make([]uint, k)

	hash := fnv.New64a()
	hash.Write(data)

	h0 := uint64(uint32(hash.Sum64()))
	h1 := uint64(uint32(hash.Sum64() >> 32))

	for i := uint32(0); i < k; i++ {
		ii := uint64(i)
		locs[i] = uint((h0 + ii*h1) % uint64(m))
	}

	return locs
}

// However, Murmur3 is 128-bit, so we can use the hash once, split technique
// and keep 64-bit addressing.
func murmurLocations(data []byte, k, m uint32) []uint {
	v1, v2 := murmur3.Sum128(data)

	positions := make([]uint, k)
	for i := range positions {
		positions[i] = uint((v1 + uint64(i)*v2) % uint64(m))
	}
	return positions
}


// The test runs many rounds, and the filter's m and k are easily adjustable.
// It builds a set of raw input data that it can re-use, which is faster but 
// uses more memory.
func TestBloomLocationUniformity(t *testing.T) {

	var m, k, rounds uint32

	m = 8
	k = 3

	rounds = 15000000

	elements := make([][]byte, rounds)

	for x := uint32(0); x < rounds; x++ {
		ctrlist := make([]uint8, 4)
		ctrlist[0] = uint8(x)
		ctrlist[1] = uint8(x >> 8)
		ctrlist[2] = uint8(x >> 16)
		ctrlist[3] = uint8(x >> 24)
		data := []byte(ctrlist)
		elements[x] = data
	}
	// fmt.Println(elements)

	fmt.Println("Willf Bloom w/4 hashes")
	chiTestBloom(m, k, rounds, elements, bloomLocations)
	fmt.Println("")

	fmt.Println("Murmur3 128-bit split to 2 64-bit hashes")
	chiTestBloom(m, k, rounds, elements, murmurLocations)
	fmt.Println("")

	fmt.Println("Willf Bloom w/2 hashes")
	chiTestBloom(m, k, rounds, elements, wfnvLocations)
	fmt.Println("")

	fmt.Println("FNV-1 64-bit split to 2 32-bit hashes")
	chiTestBloom(m, k, rounds, elements, fnvLocations)
	fmt.Println("")

	fmt.Println("FNV-1a 64-bit split to 2 32-bit hashes")
	chiTestBloom(m, k, rounds, elements, fnvLocations)
	fmt.Println("")

}

/*

This helper function ranges over the input data, applying the hashing
method fLoc(), which returns the bit locations to set in the filter.
For each location, increment a counter for that bit address.

If the Bloom Filter's location() method distributes locations uniformly
at random, a property it should inherit from its hash function, then
each bit location in the filter should end up with roughly the same
number of hits.  Importantly, the value of k should not matter.

Once the results are collected, we can run a chi squared goodness of fit
test, comparing the result histogram with the uniform distribition.
This yields a test statistic with degrees-of-freedom of m-1.

Should we desire, we can integrate the CDF of the chi square distribution
using the test statistic and the DF=m-1 to compute the p-value, i.e.,
the probability that the results were drawn from the uniform distribution.
Turns out that isn't necessary since some of the results produce a test
statistic so large.

*/
func chiTestBloom(m, k, rounds uint32, elements [][]byte, fLoc locations) {

	results := make([]uint, m)
	chi := make([]float64, m)

	var chi_statistic float64

	for _, data := range elements {
		for _, loc := range fLoc(data, k, m) {
			results[loc]++
		}
	}

	// Each element of results should contain the same value: k * rounds / m.
	// Let's run a chi-square goodness of fit and see how it fares.
	e := float64(k * rounds) / float64(m)
	for i := uint32(0); i < m; i++ {
		chi[i] = math.Pow(float64(results[i]) - e, 2.0) / e
		chi_statistic += chi[i]
	}

	fmt.Println(results)
	fmt.Println(chi)
	fmt.Println(chi_statistic)
}

func BenchmarkBloomLocations(b *testing.B) {
	for n:= 0; n < b.N; n++ {
		bloomLocations([]byte("test"), 3, 64)
	}
}
func BenchmarkMurmurLocations(b *testing.B) {
	for n:= 0; n < b.N; n++ {
		murmurLocations([]byte("test"), 3, 64)
	}
}
func BenchmarkWfnvLocations(b *testing.B) {
	for n:= 0; n < b.N; n++ {
		wfnvLocations([]byte("test"), 3, 64)
	}
}
func BenchmarkFnvLocations(b *testing.B) {
	for n:= 0; n < b.N; n++ {
		fnvLocations([]byte("test"), 3, 64)
	}
}
func BenchmarkFnvaLocations(b *testing.B) {
	for n:= 0; n < b.N; n++ {
		fnvaLocations([]byte("test"), 3, 64)
	}
}
