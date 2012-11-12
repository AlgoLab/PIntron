/*

Copyright (C) 2003-2006 Christian Kreibich <christian@whoop.org>.
Copyright (C) 2010 Yuri Pirola <yuri.pirola (*) gmail.com>.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies of the Software and its documentation and acknowledgment shall be
given in the documentation and software packages that this Software was
used.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/
#ifndef __lst_string_h
#define __lst_string_h

#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <sys/queue.h>

#define LST_STRING_END   UINT_MAX
#define LST_EMPTY_STRING UINT_MAX

typedef struct lst_string            LST_String;
typedef struct lst_stringindex       LST_StringIndex;
typedef struct lst_stringset         LST_StringSet;

/**
 * LST_StringCB - generic string callback signature.
 * @string: string passed in.
 * @data: arbitrary user data.
 */
typedef void (*LST_StringCB) (LST_String *string, void *data);


struct lst_string
{
  /* Every string has a unique id -- when appropriate
   * we don't compare string items for equality checks but
   * only look at string ids.
   */
  int                         id;

  /* Strings can be members of a set: */
  LIST_ENTRY(lst_string)      set;
  
  void                       *data;
  u_char                      data_external;

  u_int                       num_items;
  u_int                       item_size;

};


/* To implement edge-label compression, each edge is associated
 * with a string index structure, describing a substring of
 * another string by giving the start and end index (latter is
 * inclusive). Indices start at 0. For example, in string
 * "example", index (2,4) --> "amp".
 */
struct lst_stringindex
{
  /* Which string we're indexing here */
  LST_String  *string;
  
  /* Start and end indices in the string */
  u_int        start_index;

  /* The end index is kept as a pointer so
   * that it can be adjusted globally, if
   * necessary.
   */
  u_int       *end_index;
  u_int        end_index_local;

};


struct lst_stringset
{
  LIST_HEAD(slist, lst_string) members;  
  int size;
};


/**
 * lst_string_init - initializes existing string object.
 * @string: string object to initialize.
 * @data: string data to initialize with.
 * @item_size: size of a single string item, in bytes.
 * @num_items: length of string.
 *
 * The function initializes an existing string object, making it
 * use the passed data directly without copying it. It is thus faster
 * than lst_string_new, especially for e.g. tight loops.
 */
void             lst_string_init(LST_String *string, void *data,
				 u_int item_size, u_int num_items);


/**
 * lst_string_free - string destructor.
 * @string: string to clean up.
 *
 * The function cleans up all of the memory occupied by @string.
 * It is safe to call this on strings initialized using lst_string_init()
 * which are using external data, as in that case the string data
 * itself is not touched.
 */
void             lst_string_free(LST_String *string);


/**
 * lst_string_free_keep_data - string destructor, not touching string data.
 * @string: string to clean up.
 *
 * The function cleans up the memory occupied by @string without
 * releasing the actual string data, which it returns to the caller.
 *
 * Returns: actual string data.
 */
void            *lst_string_free_keep_data(LST_String *string);


/**
 * lst_string_get_length - returns number of items in string.
 * @string: string to return length of.
 * 
 * The function returns the number of items in the string. Always
 * use this function and never access any of the string members directly
 * to obtain that value.
 *
 * Returns: length of @string.
 */
u_int            lst_string_get_length(LST_String *string);


/**
 * lst_string_get_item - returns a pointer to the nth string item.
 * @string: string to look up item in.
 * @index: number of element in string to find, starting at 0.
 *
 * Returns: a pointer to the element in the string at position @index,
 * or %NULL when the index is invalid.
 */
void            *lst_string_get_item(LST_String *string, u_int index);


/**
 * lst_string_print - returns an ASCII string representation of a string.
 * @string: string to print.
 *
 * The creates an ASCII string verions of @string and returns it as
 * a pointer to static memory. The way this mapping is implemented depends
 * on the string class active for this string, see @lst_string_set_class().
 *
 * Returns: pointer to static string buffer.
 */
const char      *lst_string_print(LST_String *string);


/**
 * lst_string_item_copy - copies an item from one string to another.
 * @src: source string.
 * @src_index: item in source string, starting at 0.
 * @dst: destination string.
 * @dst_index: item to copy to.
 *
 * The function copies the item found at @src_index in @src into the
 * item @dst_index of string @dst.
 */
void             lst_string_item_copy(LST_String *src, u_int src_index,
				      LST_String *dst, u_int dst_index);


/**
 * lst_string_eq - string item comparison.
 * @s1: first string.
 * @item1: item in @s1.
 * @s2: second string.
 * @item2: item in @s2.
 *
 * The function compares the items specified via the input parameters.
 * The way this is implemented depends on the string class for the
 * strings involved, see lst_string_set_class().
 *
 * Returns: value > 0 when equal, 0 otherwise.
 */
static inline
int              lst_string_eq(LST_String *s1, u_int item1,
			       LST_String *s2, u_int item2);

static inline
int
lst_string_eq(LST_String *s1, u_int item1,
				  LST_String *s2, u_int item2)
{
  if (!s1 || !s2 || item1 >= s1->num_items || item2 >= s2->num_items)
    return 0;

  /* Treat the end-of-string markers separately: */
  if (item1 == s1->num_items - 1 || item2 == s2->num_items - 1) {
	 if (item1 == s1->num_items - 1 && item2 == s2->num_items - 1) {
		if (s1 == s2) {
		  return 1;
		} else {
		  return 0;
		}
	 } else {
		return 0;
	 }
  }

  return *((char*)lst_string_get_item(s1, item1)) == *((char*)lst_string_get_item(s2, item2));
}


/**
 * lst_string_items_common - find string overlap at specific indices.
 * @s1: first string.
 * @off1: item in @s1.
 * @s2: second string.
 * @off2: item in @s2.
 * @max_len: maximum number of items to compare.
 *
 * The function compares items in @s1 and @s2 from the given offsets,
 * counting how many are equal. The way the comparison works depends on
 * the string class active for the strings involved, see lst_string_set_class().
 *
 * Returns: number of identical items.
 */
u_int            lst_string_items_common(LST_String *s1, u_int off1,
					 LST_String *s2, u_int off2,
					 u_int max_len);


/**
 * lst_string_print_hex - string printer printing string as hex byte data.
 * @index: string region to print.
 * 
 * For convenience, this is a printer implementation that returns a hex
 * representation of the data contained in the string. You can call
 * this function up to three times from within the same printf()-like
 * function.
 *
 * Returns: hex representation in statically allocated string, copy
 * this if you need to keep it around.
 */
char *           lst_string_print_hex(LST_StringIndex *index);



/**
 * lst_string_index_init - initializes a string index.
 * @index: index initialized.
 *
 * The function initializes a string index. Used internally.
 */
void             lst_string_index_init(LST_StringIndex *index);

/**
 * lst_string_index_copy - copies a string index.
 * @src: source index.
 * @dst: destination index.
 *
 * Used internally.
 */
void             lst_string_index_copy(LST_StringIndex *src, LST_StringIndex *dst);


/**
 * lst_stringset_new - creates a new stringset.
 * 
 * The function creates a new stringset. Stringsets are the way you
 * pass multiple strings to an algorithm. You basically create a string
 * set, then add strings to the set, pass the set to an algorithm, and
 * clean up the set afterwards.
 *
 * Returns: new, empty set.
 */
LST_StringSet   *lst_stringset_new(void);


/**
 * lst_stringset_add - adds a string to a string set.
 * @set: set to add string to.
 * @string: string to add.
 *
 * The function adds a string to the set.
 */
void             lst_stringset_add(LST_StringSet *set, LST_String *string);


/**
 * lst_stringset_remove - removes a string from a string set.
 * @set: set to remove string from.
 * @string: string to remove.
 *
 * The function removes a string from a string set.
 */
void             lst_stringset_remove(LST_StringSet *set, LST_String *string);


/**
 * lst_stringset_foreach - calls a callback for each string in the set.
 * @set: set to iterate.
 * @callback: callback to call for each item.
 * @user_data: arbitrary data passed through to @callback.
 *
 * The function calls @callback for each string in @set, passing it
 * that string and @user_data.
 */
void             lst_stringset_foreach(LST_StringSet *set, LST_StringCB callback, void *user_data);


/**
 * lst_stringset_free - cleans up a string set.
 * @set: set to clean up.
 *
 * The function releases all the memory claimed by @set, including
 * all the strings it contains.
 */
void             lst_stringset_free(LST_StringSet *set);

#endif
