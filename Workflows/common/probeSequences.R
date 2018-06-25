wl.file = function(lines, o.file)
{
   f.conn = o.file %>% file
   open(f.conn, open="w")
   if (isOpen(f.conn))
   {   
      writeLines(lines, f.conn)
      close(f.conn)
      cat("Closing:", o.file, "\n")
      return(T)
   } else warning("Connection to ", restFile, " could not be established")

   return(F)
}
