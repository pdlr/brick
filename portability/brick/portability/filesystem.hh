/**
***************************************************************************
* @file brick/portability/filesystem.hh
*
* Header file declaring portability routines for dealing with filesystems.
*
* Copyright (C) 2003-2011, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/


#ifndef BRICK_PORTABILITY_FILESYSTEM_HH
#define BRICK_PORTABILITY_FILESYSTEM_HH

#include <string>
#include <vector>

// Anonymous namespace for stuff local to this file.
namespace brick {

  namespace portability {

    /** 
     * This function returns the preferred path delimiter for the
     * operating system.  For example, on Linux operating systems, it
     * will return "/".
     * 
     * @return The return value is a string containing the preferred
     * path delimeter.
     */
    const std::string& pathDelimiter();


    /** 
     * This function returns a vector of strings containing all of the
     * recognized path delimiters for the operating system.  The first
     * element of this vector will always be equal to the string by
     * function pathDelimiter().  Some operating systems may use more
     * than one path delimeter, and in this case, the vector will have
     * more than one element.  For example, on Windows, the returned
     * vector has two elements: "\", and "/".
     * 
     * @return The return value is a vector of recognized path
     * delimiters.
     */
    const std::vector<std::string>& pathDelimiters();


    /** 
     * This function returns the string that separates a filename
     * extension from the rest of the filename.  On all operating
     * systems so far, this delimiter is ".".  For example, in the
     * filename /foo/bar/baz.dat, the filename ("baz") is separated
     * from the extension ("dat") by a single ".".
     * 
     * @return The return value is a string containing the preferred 
     * delimiter.
     */
    const std::string& extensionDelimiter();


    /** 
     * This function returns a bool indicating whether or not the
     * specified path is a directory.  If the path is a symbolic link,
     * the return value is currently unspecified, but will eventually be
     * true iff the link points to a directory.
     * 
     * @param path This argument is the filename to evaluate.
     * 
     * @return The return value is true if path refers to a directory,
     * false otherwise.
     */
    bool
    isDirectory(const std::string& path);
    

    /** 
     * This function combines parts of a path using the appropriate
     * path delimiter.  Here are some examples:
     *
     * @code
     *   // These calls return "/foo/bar/baz" 
     *   joinPath("/foo", "bar/baz");
     *   joinPath("/foo/", "bar/baz");
     *
     *   // This call calls returns "/bar/baz"
     *   joinPath("/foo", "/bar/baz");
     * @endcode
     *
     * @param part0 This argument is the first part of the path.
     * 
     * @param part1 This argument is the second part of the path.
     * 
     * @return The return value is the two parts, joined together 
     * in a way that is appropriate for the host operating system.
     */
    std::string
    joinPath(const std::string& part0, const std::string& part1);


    /** 
     * This function lists the contents of a directory.  The returned
     * vector is not sorted.
     * 
     * @param directoryName This argument is the directory to be
     * examined.
     * 
     * @param fullPath This argument specifies whether the returned
     * filenames and directory names should be prepended with the name
     * of the directory that was listed.
     * 
     * @return The return value is a vector of filenames and directory
     * names.
     */
    std::vector<std::string>
    listDirectory(const std::string& directoryName, bool fullPath);


    /** 
     * This function accepts a path returns a pair of strings in which
     * the first element is the directory name and the second is the
     * filename.  For example:
     *
     *   splitPath("/foo/bar.baz")
     *
     * returns
     *
     *   {"/foo/", "bar.baz"}
     *
     * Also,
     *
     *   splitPath("bar.baz")
     *
     * returns
     *
     *   {"", "bar.baz"}
     *
     *  and
     * 
     *   splitPath("/foo/")
     *
     * returns
     *
     *   {"/foo/", ""}
     *
     * @return The return value is a pair containing first the directory
     * name, and second the file name.
     */
    std::pair<std::string, std::string>
    splitPath(const std::string& path);
    
  } // namespace portability
  
} // namespace brick

#endif /* ifndef BRICK_PORTABILITY_FILESYSTEM_HH */
