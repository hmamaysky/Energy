require(data.table)

#' Some of the text munging.
DriverEneryData <- function() {

  ## the text summaries
  ws <- getTopicWordLists()

  ssent <- getSampleSentences()

}

#' Sample sentences.
getSampleSentences <- function() {

  dname <- '~/Dropbox/energy drivers/data/textual measures/'
  wname <- '~/GoogleDrive/writing/WP/energy'

  ## read the sample articles
  ssent <- fread(file.path(dname,'info_sample_oil_articles.csv'))

  setnames(ssent,c('sentiment','article_topic','headline'),c('Sentiment','Topic','Headline'))

  ## change topics names
  ssent[, Topic := gsub('Topic [0-9]( )*:( )','',Topic)]

  ## get the first and last 2 articles
  samplesent <- ssent[,.SD[c(1,2,nrow(.SD)-c(0,1))],by=Topic][,.(Topic,Sentiment,Headline)]

  write.csv(samplesent,file=file.path(wname,paste0('sample-sentences-',Sys.Date(),'.csv')))

  samplesent
}


#' The topic word lists.
getTopicWordLists <- function() {

  dname <- '~/Dropbox/energy drivers/data/textual measures/'
  wname <- '~/GoogleDrive/writing/WP/energy'

  require(XLConnect)
  wb <- loadWorkbook(file.path(dname,'oil_words_clustering_C.xlsx'))
  ws <- readWorksheet(wb,sheet='long form',header=TRUE) %>% as.data.table

  ## output the word lists for each topic
  wordlists <- ws[order(Total,decreasing=TRUE),
     .(WordList=
         paste0(sapply(1:20,function(ii) sprintf('%s (%s)',Word[[ii]],format(Total[[ii]],big.mark=','))),
                collapse=', ')
       ), by=Topic]
  write.csv(wordlists, file=file.path(wname,paste0('wordlists-',Sys.Date(),'.csv')))

  wordlists
}
