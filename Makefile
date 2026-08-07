# ==============================================================================
# Removal_Study_BLIMMP pipeline Makefile
# ==============================================================================
# Each pipeline stage is a Make target. Make decides what needs to (re)run by
# comparing timestamps of ".done" stamp files -- not the many individual output
# files each stage produces, since there are too many of those to list cleanly.
#
# Usage:
#   make                # run the whole pipeline in order, skipping finished stages
#   make subsample      # run only up through the subsample stage
#   make clean-hmmer    # remove stamps + outputs for the hmmer stage (forces rerun)
#   make status         # show which stages are done / pending
#
# How targets get "done": each stage submits a SLURM job (via sbatch) and BLOCKS
# until it finishes (see `run_and_wait` below), then touches its .done stamp.
# This keeps `make` itself simple -- run it in a login shell, a screen/tmux
# session, or as its own lightweight sbatch job if you want it to survive
# logout.
# ==============================================================================

SHELL := /bin/bash
BASE_DIR := $(shell pwd)
STAMP_DIR := $(BASE_DIR)/.make_stamps

.PHONY: all status clean clean-tantan clean-prodigal clean-busco clean-subsample \
        clean-hmmer clean-dram clean-dram-distill clean-anvio clean-metabolic

all: busco hmmer dram-distill anvio-estimate metabolic

$(STAMP_DIR):
	mkdir -p $(STAMP_DIR) logs

# ------------------------------------------------------------------------------
# Helper: submit an sbatch script and block until it completes.
# Fails (non-zero exit) if the job ends in a non-COMPLETED state, so Make will
# correctly stop the chain rather than touching a .done stamp for a failed job.
# ------------------------------------------------------------------------------
define run_and_wait
	jobid=$$(sbatch --parsable $(1)); \
	echo "Submitted $(1) as job $$jobid, waiting..."; \
	while true; do \
		state=$$(sacct -j $$jobid.0 --format=State --noheader 2>/dev/null | head -1 | tr -d ' '); \
		if [[ -z "$$state" ]]; then sleep 15; continue; fi; \
		case "$$state" in \
			COMPLETED) echo "Job $$jobid COMPLETED"; break ;; \
			FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL) \
				echo "Job $$jobid ended in state $$state -- aborting." >&2; exit 1 ;; \
			*) sleep 30 ;; \
		esac; \
	done
endef

# ------------------------------------------------------------------------------
# Stage 1: tantan masking + prodigal gene calling
# ------------------------------------------------------------------------------
$(STAMP_DIR)/tantan_prodigal.done: run_tantan_prodigal.sh | $(STAMP_DIR)
	$(call run_and_wait,run_tantan_prodigal.sh)
	touch $@

tantan-prodigal: $(STAMP_DIR)/tantan_prodigal.done

# ------------------------------------------------------------------------------
# Stage 2: BUSCO completeness (depends on nothing upstream in this chain --
# runs against raw/masked genomes directly)
# ------------------------------------------------------------------------------
$(STAMP_DIR)/busco.done: run_busco.sh $(STAMP_DIR)/tantan_prodigal.done | $(STAMP_DIR)
	$(call run_and_wait,run_busco.sh)
	touch $@

busco: $(STAMP_DIR)/busco.done

# ------------------------------------------------------------------------------
# Stage 3: ORF subsampling (needs prodigal's *_ORFs.faa from stage 1)
# ------------------------------------------------------------------------------
$(STAMP_DIR)/subsample.done: run_subsample_orfs.sh subsample_orfs.py $(STAMP_DIR)/tantan_prodigal.done | $(STAMP_DIR)
	$(call run_and_wait,run_subsample_orfs.sh)
	touch $@

subsample: $(STAMP_DIR)/subsample.done

# ------------------------------------------------------------------------------
# Stage 4: full-genome HMMER (split -> array search -> merge)
# ------------------------------------------------------------------------------
$(STAMP_DIR)/hmm_split.done: split_hmm_profiles.sh $(STAMP_DIR)/tantan_prodigal.done | $(STAMP_DIR)
	$(call run_and_wait,split_hmm_profiles.sh)
	touch $@

$(STAMP_DIR)/hmmsearch.done: run_hmmsearch_array.sh $(STAMP_DIR)/hmm_split.done | $(STAMP_DIR)
	$(call run_and_wait,run_hmmsearch_array.sh)
	touch $@

$(STAMP_DIR)/hmmer_merge.done: merge_hmmsearch_results.sh $(STAMP_DIR)/hmmsearch.done | $(STAMP_DIR)
	$(call run_and_wait,merge_hmmsearch_results.sh)
	touch $@

hmmer: $(STAMP_DIR)/hmmer_merge.done

# ------------------------------------------------------------------------------
# Stage 5: DRAM annotate_genes on every subsample
# ------------------------------------------------------------------------------
$(STAMP_DIR)/dram_annotate.done: run_dram_array.sh $(STAMP_DIR)/subsample.done | $(STAMP_DIR)
	$(call run_and_wait,run_dram_array.sh)
	touch $@

dram-annotate: $(STAMP_DIR)/dram_annotate.done

# ------------------------------------------------------------------------------
# Stage 6: DRAM distill (needs annotate_genes output)
# ------------------------------------------------------------------------------
$(STAMP_DIR)/dram_distill.done: run_dram_distill_array.sh $(STAMP_DIR)/dram_annotate.done | $(STAMP_DIR)
	$(call run_and_wait,run_dram_distill_array.sh)
	touch $@

dram-distill: $(STAMP_DIR)/dram_distill.done

# ------------------------------------------------------------------------------
# Stage 7: anvi'o contigs-db + KOfam build (needs subsamples)
# ------------------------------------------------------------------------------
$(STAMP_DIR)/anvio_build.done: run_anvio_build_array.sh make_external_gene_calls.py $(STAMP_DIR)/subsample.done | $(STAMP_DIR)
	$(call run_and_wait,run_anvio_build_array.sh)
	touch $@

anvio-build: $(STAMP_DIR)/anvio_build.done

# ------------------------------------------------------------------------------
# Stage 8: anvi'o metabolism estimation (needs the built contigs.db)
# ------------------------------------------------------------------------------
$(STAMP_DIR)/anvio_estimate.done: run_anvio_estimate_array.sh $(STAMP_DIR)/anvio_build.done | $(STAMP_DIR)
	$(call run_and_wait,run_anvio_estimate_array.sh)
	touch $@

anvio-estimate: $(STAMP_DIR)/anvio_estimate.done

# ------------------------------------------------------------------------------
# Stage 9: METABOLIC on every subsample (independent of DRAM/anvi'o -- only
# needs the subsampled ORF faa files)
# ------------------------------------------------------------------------------
$(STAMP_DIR)/metabolic.done: run_metabolic_array.sh $(STAMP_DIR)/subsample.done | $(STAMP_DIR)
	$(call run_and_wait,run_metabolic_array.sh)
	touch $@

metabolic: $(STAMP_DIR)/metabolic.done

# ------------------------------------------------------------------------------
# Status / cleaning
# ------------------------------------------------------------------------------
status:
	@echo "Pipeline status:"
	@for stage in tantan_prodigal busco subsample hmm_split hmmsearch hmmer_merge dram_annotate dram_distill anvio_build anvio_estimate metabolic; do \
		if [[ -f "$(STAMP_DIR)/$$stage.done" ]]; then \
			echo "  [x] $$stage  (done $$(date -r $(STAMP_DIR)/$$stage.done '+%Y-%m-%d %H:%M'))"; \
		else \
			echo "  [ ] $$stage"; \
		fi; \
	done

clean-tantan:
	rm -f $(STAMP_DIR)/tantan_prodigal.done

clean-busco:
	rm -f $(STAMP_DIR)/busco.done

clean-subsample:
	rm -f $(STAMP_DIR)/subsample.done

clean-hmmer:
	rm -f $(STAMP_DIR)/hmm_split.done $(STAMP_DIR)/hmmsearch.done $(STAMP_DIR)/hmmer_merge.done

clean-dram:
	rm -f $(STAMP_DIR)/dram_annotate.done

clean-dram-distill:
	rm -f $(STAMP_DIR)/dram_distill.done

clean-anvio:
	rm -f $(STAMP_DIR)/anvio_build.done $(STAMP_DIR)/anvio_estimate.done

clean-metabolic:
	rm -f $(STAMP_DIR)/metabolic.done

clean:
	rm -rf $(STAMP_DIR)
