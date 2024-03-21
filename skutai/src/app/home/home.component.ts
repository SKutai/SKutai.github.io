import { CommonModule } from '@angular/common';
import { AfterViewInit, Component, OnDestroy, OnInit } from '@angular/core';
import { Observable, Subscription, interval } from 'rxjs';

@Component({
  selector: 'app-home',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.scss'],
})
export class HomeComponent implements OnInit, AfterViewInit, OnDestroy {
  private subscription1: Subscription | undefined;
  private subscription2: Subscription | undefined;

  interval1$: Observable<number> = interval(75);
  interval2$: Observable<number> = interval(1000);

  fullGreeting: string = "Hi there! |Thanks for visiting my website!\n| Feel free to use any of the buttons in the top right to look around.\n| Or just watch the city go by for a while |- that's cool too.";
  cursor: string = "█";
  greeting: string = "█";
  typing: boolean = true;
  characterIndex: number = 0;
  skipCount: number = 0;
  startingDelay: number = 1000;
  endingDelay: number = 10000;
  doneGreeting: boolean = false;

  constructor() {

  }

  ngOnInit(): void {

  }

  ngAfterViewInit(): void {
    // wait for a moment then type out the greeting
    setTimeout(() => {
      this.subscription1 = this.interval1$.subscribe((x) => {
        if (this.skipCount <= 0) {
          if (this.characterIndex <= this.fullGreeting.length) {
            let takeExtraChar = false;
            let PutCursorOnNewLine = false;

            if (this.fullGreeting[this.characterIndex] == "\n") {
              PutCursorOnNewLine = true;
            }

            if (this.fullGreeting[this.characterIndex - 1] == "\n") {
              // insert a space in the full greeting
              this.fullGreeting = this.fullGreeting.slice(0, this.characterIndex - 1) + " " + this.fullGreeting.slice(this.characterIndex - 1);
              this.characterIndex++;
              takeExtraChar = true;
            }

            if (takeExtraChar) {
              this.greeting = this.fullGreeting.slice(0, this.characterIndex + 1) + this.cursor;
            }
            else if (PutCursorOnNewLine) {
              this.greeting = this.fullGreeting.slice(0, this.characterIndex) + " \n " + this.cursor;
            }
            else {
              this.greeting = this.fullGreeting.slice(0, this.characterIndex) + this.cursor;
            }

            // if the next character is |
            if (this.fullGreeting[this.characterIndex + 1] == "|") {
              // remove the |
              this.fullGreeting = this.fullGreeting.replace("|", "");
              this.pauseTyping(10);
            }

            this.characterIndex++;
          }
          else {
            this.typing = false;
            // wait a moment before declaring the greeting done
            setTimeout(() => {
              this.doneGreeting = true;
            }, this.endingDelay);
          }
        }
        else {
          this.skipCount--;
        }
      });

      this.subscription2 = this.interval2$.subscribe((x) => {
        if (!this.typing) {
          if (x % 2 == 0 && !this.greeting.includes(this.cursor)) {
            this.greeting = this.fullGreeting.trim() + this.cursor;
          }
          else {
            this.greeting = this.greeting.replaceAll("█", " ");
          }
        }
      });


    }, this.startingDelay);
  }

  pauseTyping(skipCount: number) {
    this.skipCount = skipCount;
  }

  ngOnDestroy(): void {
    if (this.subscription1) {
      this.subscription1.unsubscribe();
    }
    if (this.subscription2) {
      this.subscription2.unsubscribe();
    }
  }

}
