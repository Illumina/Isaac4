################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2014 Illumina, Inc.
## All rights reserved.
##
## This software is provided under the terms and conditions of the
## GNU GENERAL PUBLIC LICENSE Version 3
##
## You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
## along with this program. If not, see
## <https://github.com/illumina/licenses/>.
##
################################################################################
##
## file UnpackReferenceSubmake.mk
##
## brief Defines appropriate rules
##
## author Roman Petrovski
##
################################################################################

MAKEFILES_DIR:=@iSAAC_HOME@@iSAAC_FULL_DATADIR@/makefiles

# Import the global configuration
include $(MAKEFILES_DIR)/common/Config.mk

include $(MAKEFILES_DIR)/common/Sentinel.mk

# Import the logging functionalities
include $(MAKEFILES_DIR)/common/Log.mk

# Import the debug functionalities
include $(MAKEFILES_DIR)/common/Debug.mk

include $(MAKEFILES_DIR)/reference/Config.mk

ifeq (yes,$(MOVABLE))
ROOT_PATH:=.
else
ROOT_PATH:=$(CURDIR)
endif


UNPACKED_SORTED_REFERENCE_XML:=$(TEMP_DIR)/sorted-reference.xml
ORIGINAL_GENOME_FASTA:=$(notdir $(shell $(XSLTPROC) $(GET_ANY_FASTA_PATH_XSL) $(UNPACKED_SORTED_REFERENCE_XML)))

UNPACKED_GENOME_FILE:=$(TEMP_DIR)/$(ORIGINAL_GENOME_FASTA)
GENOME_FILE:=$(ROOT_PATH)/$(ORIGINAL_GENOME_FASTA)

SORTED_REFERENCE_XML:=sorted-reference.xml

get_format_version=$(shell $(XSLTPROC) $(GET_FORMAT_VERSION_XSL) $(UNPACKED_SORTED_REFERENCE_XML))

ifneq ($(CURRENT_REFERENCE_FORMAT_VERSION),$(get_format_version))
$(error "Unsupported packed reference format $(get_format_version). Version $(CURRENT_REFERENCE_FORMAT_VERSION) is required")
endif

$(GENOME_FILE): $(UNPACKED_GENOME_FILE)
	$(CMDPREFIX) $(CP) $< $(SAFEPIPETARGET)

$(CONTIGS_XML): $(GENOME_FILE) $(TEMP_DIR)/.sentinel $(UNPACKED_SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(PRINT_CONTIGS) -g $(GENOME_FILE) --original-metadata $(UNPACKED_SORTED_REFERENCE_XML) >$(SAFEPIPETARGET)

$(SORTED_REFERENCE_XML): $(CONTIGS_XML)
	$(CMDPREFIX) $(MERGE_REFERENCES) --make-absolute-paths no $(foreach part, $^, -i '$(part)') -o $(SAFEPIPETARGET)

