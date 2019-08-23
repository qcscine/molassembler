/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Base 64 encoding and decoding between strings and vectors of uint_8
 */

#ifndef INCLUDE_BASE_64_ENCODING_H
#define INCLUDE_BASE_64_ENCODING_H

#include <string>
#include <vector>
#include <stdexcept>

namespace Scine {

namespace base64 {

/* Adapted from public domain licensed code at
 * https://en.wikibooks.org/wiki/Algorithm_Implementation/Miscellaneous/Base64
 */

// Lookup table
// If you want to use an alternate alphabet, change the characters here
constexpr char encodeLookup[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
constexpr char padCharacter = '=';

/** @brief Encode binary data as a string
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param inputBuffer binary data
 */
std::string encode(const std::vector<std::uint8_t>& inputBuffer) {
	std::string encodedString;
	encodedString.reserve(
    (
      (inputBuffer.size() / 3)
      + (inputBuffer.size() % 3 > 0)
    ) * 4
  );
	long temp;
	std::vector<std::uint8_t>::const_iterator cursor = inputBuffer.begin();

	for(size_t idx = 0; idx < inputBuffer.size()/3; idx++) {
		temp  = (*cursor++) << 16; //Convert to big endian
		temp += (*cursor++) << 8;
		temp += (*cursor++);
		encodedString.append(1, encodeLookup[(temp & 0x00FC0000) >> 18]);
		encodedString.append(1, encodeLookup[(temp & 0x0003F000) >> 12]);
		encodedString.append(1, encodeLookup[(temp & 0x00000FC0) >> 6 ]);
		encodedString.append(1, encodeLookup[(temp & 0x0000003F)      ]);
	}

	switch(inputBuffer.size() % 3) {
	case 1:
		temp  = (*cursor++) << 16; //Convert to big endian
		encodedString.append(1, encodeLookup[(temp & 0x00FC0000) >> 18]);
		encodedString.append(1, encodeLookup[(temp & 0x0003F000) >> 12]);
		encodedString.append(2, padCharacter);
		break;
	case 2:
		temp  = (*cursor++) << 16; //Convert to big endian
		temp += (*cursor++) << 8;
		encodedString.append(1, encodeLookup[(temp & 0x00FC0000) >> 18]);
		encodedString.append(1, encodeLookup[(temp & 0x0003F000) >> 12]);
		encodedString.append(1, encodeLookup[(temp & 0x00000FC0) >> 6 ]);
		encodedString.append(1, padCharacter);
		break;
	}

	return encodedString;
}

/** @brief Decode base 64 string data to binary
 *
 * @complexity{@math{\Theta(N)}}
 *
 * @param input base 64 string
 */
std::vector<std::uint8_t> decode(const std::string& input) {
	if(input.length() % 4) {
		throw std::runtime_error("Invalid base64 encoding!");
  }

	size_t padding = 0;

	if(input.length()) {
		if(input[input.length()-1] == padCharacter) {
			++padding;
    }

		if(input[input.length()-2] == padCharacter) {
			++padding;
    }
	}

	//Setup a vector to hold the result
	std::vector<std::uint8_t> decodedBytes;
	decodedBytes.reserve(((input.length()/4)*3) - padding);
	long temp=0; //Holds decoded quanta
	std::string::const_iterator cursor = input.begin();

	while(cursor < input.end()) {
		for(size_t quantumPosition = 0; quantumPosition < 4; quantumPosition++) {
			temp <<= 6;
			if(*cursor >= 0x41 && *cursor <= 0x5A) // This area will need tweaking if
				temp |= *cursor - 0x41;		              // you are using an alternate alphabet
			else if(*cursor >= 0x61 && *cursor <= 0x7A)
				temp |= *cursor - 0x47;
			else if(*cursor >= 0x30 && *cursor <= 0x39)
				temp |= *cursor + 0x04;
			else if(*cursor == 0x2B)
				temp |= 0x3E; //change to 0x2D for URL alphabet
			else if(*cursor == 0x2F)
				temp |= 0x3F; //change to 0x5F for URL alphabet
			else if(*cursor == padCharacter) {
        // Pad the rest
				switch(input.end() - cursor) {
				case 1: //One pad character
					decodedBytes.push_back((temp >> 16) & 0x000000FF);
					decodedBytes.push_back((temp >> 8 ) & 0x000000FF);
					return decodedBytes;
				case 2: //Two pad characters
					decodedBytes.push_back((temp >> 10) & 0x000000FF);
					return decodedBytes;
				default:
					throw std::runtime_error("Invalid Padding in Base 64!");
				}
			} else {
				throw std::runtime_error("Non-Valid Character in Base 64!");
      }

			++cursor;
		}

		decodedBytes.push_back((temp >> 16) & 0x000000FF);
		decodedBytes.push_back((temp >> 8 ) & 0x000000FF);
		decodedBytes.push_back((temp      ) & 0x000000FF);
	}

	return decodedBytes;
}

} // namespace base64

} // namespace Scine

#endif
