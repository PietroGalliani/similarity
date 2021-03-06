-- MySQL dump 10.13  Distrib 5.1.50, for unknown-linux-gnu (x86_64)
--
-- Host: localhost    Database: go
-- ------------------------------------------------------
-- Server version	5.1.50-community

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `graph_path`
--

DROP TABLE IF EXISTS `graph_path`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `graph_path` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `term1_id` int(11) NOT NULL,
  `term2_id` int(11) NOT NULL,
  `relationship_type_id` int(11) DEFAULT NULL,
  `distance` int(11) DEFAULT NULL,
  `relation_distance` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `graph_path0` (`id`),
  KEY `relationship_type_id` (`relationship_type_id`),
  KEY `graph_path1` (`term1_id`),
  KEY `graph_path2` (`term2_id`),
  KEY `graph_path3` (`term1_id`,`term2_id`),
  KEY `graph_path4` (`term1_id`,`distance`),
  KEY `graph_path5` (`term1_id`,`term2_id`,`relationship_type_id`),
  KEY `graph_path6` (`term1_id`,`term2_id`,`relationship_type_id`,`distance`,`relation_distance`),
  KEY `graph_path7` (`term2_id`,`relationship_type_id`),
  KEY `graph_path8` (`term1_id`,`relationship_type_id`)
) ENGINE=MyISAM AUTO_INCREMENT=1221377 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2015-08-16  1:44:01
